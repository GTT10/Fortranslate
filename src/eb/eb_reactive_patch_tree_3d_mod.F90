module eb_reactive_patch_tree_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite,ieee_next_after
  use, intrinsic :: iso_fortran_env, only: int64
  use rotating_frame_3d_mod, only: rotating_frame_3d,convert_rotating_energy_3d, &
    rotating_potential_3d,rotating_momentum_source_3d
  use eb_held_funding_queue_3d_mod, only: pack_eb_held_funding_faces_3d,physical_eb_held_funding_debt_3d, &
    pack_eb_held_funding_increment_3d,eb_held_funding_context_3d, &
    valid_eb_held_funding_context_3d,move_eb_held_funding_context_3d, &
    initialize_eb_held_funding_context_3d,copy_eb_held_funding_context_3d, &
    set_eb_held_funding_clock_faces_3d,close_eb_held_funding_clock_faces_3d, &
    read_eb_held_funding_open_faces_3d,take_eb_held_funding_target_3d,take_eb_held_funding_bucket_3d, &
    read_eb_held_funding_stock_3d,queue_eb_held_funding_receipt_3d
  use precision_mod, only: dp
  use constant_volume_reactor_mod, only: forced_reactor_integrator,advance_forced_reactor_interval
  use dual_energy_mod, only: dual_energy_policy,reset_nasa7_dual_energy_patch,reset_nasa7_dual_energy_cell, &
    initialize_nasa7_internal_energy_patch
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use simulation_config_reactive_3d_mod, only: mach_inlet_ghost_3d,valid_reactive_boundaries_3d
  use reactive_3d_mod, only: advance_reactive_chemistry_3d,react_sdc_state_3d
  use gas_transport_mod, only: gas_transport_species
  use eb_reactive_transport_3d_mod, only: advance_reactive_eb_transport_patch_euler_3d, &
    reactive_eb_isothermal_timestep_3d,compute_sdc_eb_transport_3d,compute_sdc_eb_transport_halo_3d, &
    reactive_eb_transport_patch_timestep_3d,reactive_eb_halo_bounds_3d
  use reactive_eb_cfl_3d_mod, only: reactive_eb_signal_speeds_3d,compute_reactive_eb_euler_cfl_timestep_3d
  use eb_reactive_hydro_3d_mod, only: advance_reactive_eb_euler_3d,advance_reactive_eb_plm_ctu_3d, &
    compute_reactive_eb_flattening_3d,eb_plm_boundary_3d,compute_sdc_eb_hydro_3d,limit_reactive_eb_flux_correction_3d, &
    correct_reactive_eb_spatial_hydro_flux_3d
  use eb_wall_heat_source_3d_mod, only: compute_eb_wall_transport_source_3d
  use eb_mol_hydro_3d_mod, only: compute_mol_eb_hydro_3d,eb_mol_diagnostics_3d
  use eb_reactive_wall_flux_3d_mod, only: reactive_eb_slip_wall_source_3d
  use elementary_kinetics_mod, only: elementary_reaction
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_tagging_3d_mod, only: hot_cell,build_tagged_tree_3d
  use eb_reactive_redistribution_3d_mod, only: state_redistribution_graph_slopes_3d,limit_state_recipient_slopes, &
    state_redistribution_boundary_nodes_3d,state_redist_register_3d,compute_state_redist_level_stage_3d, &
    accumulate_state_redist_register_3d,accumulate_state_redist_coarse_register_3d,assemble_state_redist_interface_rereflux_3d, &
    initialize_state_redist_tree_register_3d,accumulate_state_redist_mol_registration_stage_3d, &
    initialize_flux_rereflux_register_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d,periodic_eb_faces_match
  use eb_reactive_two_level_3d_mod, only: eb_reactive_patch_3d,matching_eb_metrics_3d
  use eb_reactive_reflux_3d_mod, only: compute_reactive_eb_reactor_bias_3d,admit_reactive_eb_reactor_energy_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d,valid_eb_integrated_flux_3d, &
    restrict_eb_integrated_flux_3d,reconcile_reactive_eb_peer_flux_3d,reflux_reactive_eb_state_3d, &
    compute_reactive_eb_reflux_increment_3d,reconcile_reactive_eb_interface_budget_3d, &
    supply_reactive_eb_interface_donors_3d,supply_reactive_eb_endpoint_donor_3d, &
    supply_reactive_eb_held_eos_donor_3d,supply_reactive_eb_gas_eos_donor_3d, &
    reactive_eb_funding_admission,compute_reactive_eb_held_projection,compute_accepted_eb_held_clock_prefixes_3d, &
    reactive_eb_exchange_interval_3d,add_eb_rereflux_to_fine_3d
  use eb_reactive_transfer_3d_mod, only: recover_reactive_eb_cell_3d, prolong_reactive_eb_state_pcm_3d, &
      restrict_reactive_eb_state_3d, &
    overlay_reactive_eb_peer_halo_3d,plan_eb_linear_halo_3d,interpolate_eb_linear_halo_3d, &
    reconstruct_eb_pressure_state_3d,interpolate_eb_pressure_time_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  use reactive_transport_3d_mod, only: valid_constant_transport_coefficients
  use eb_reactive_redistribution_3d_mod, only: post_eb_spatial_density_rate_3d
  implicit none
  private
  type, public :: eb_reactive_tree_3d
    type(amr_patch_tree_3d) :: topology
    type(eb_reactive_patch_3d), allocatable :: nodes(:)
  end type
  ! Accepted native Godunov clock data. Source evaluations and hydro passes
  ! are distinct nonlinear stages; a frozen mean cannot replace either.
  type, public :: eb_godunov_source_stage_3d
    integer :: index=-1
    logical :: coarse_registers_accumulated=.false.
    logical :: auxiliary_accepted=.false.,provider_accepted=.false.
    integer :: provider_kind=0
    real(dp) :: provider_limiter=1
    type(amr_tree_state_3d) :: provider_rate,provider_halos
    type(eb_integrated_flux_3d), allocatable :: provider_faces(:)
    type(amr_tree_state_3d) :: gas,temperature,internal_energy,external,auxiliary,rate,scalar_rate
    type(eb_integrated_flux_3d), allocatable :: faces(:),scalar_faces(:)
    type(state_redist_register_3d), allocatable :: registers(:),scalar_registers(:)
  end type eb_godunov_source_stage_3d
  ! Observed packets from a refused ordinary native MOL clock. This is not
  ! an accepted interval: future stage fields stay unallocated, and the caller
  ! retains its gas/H/ledgers. Recovery must rebuild the native spatial family.
  type, public :: eb_native_mol_source_prefix_3d
    integer :: depth=-1,iterations=0,accepted_source_stages=0
    real(dp) :: interval=0,force(3)=0,acceleration(3)=0,initial_energy_correction=0
    logical :: separate_momentum=.false.,consistent_velocity=.false.,paired_pressure=.false.,connected_directions=.false.
    logical :: chemistry=.false.,transport_flags(4)=.false.
    logical :: has_auxiliary_source=.false.,auxiliary_has_internal_energy=.false.,has_frame=.false.
    logical :: has_source_funding=.false.
    type(dual_energy_policy) :: policy
    type(rotating_frame_3d) :: frame
    type(amr_patch_tree_3d) :: topology
    type(amr_tree_state_3d) :: initial_gas,initial_temperature,initial_internal_energy,initial_reaction_source
    type(amr_tree_state_3d) :: source_funding_stock
    type(eb_geometry_3d), allocatable :: geometry(:)
    type(amr_tree_state_3d), allocatable :: external_stages(:),auxiliary_stages(:)
  end type eb_native_mol_source_prefix_3d
  type, public :: eb_godunov_hydro_stage_3d
    integer :: index=-1
    logical :: coarse_registers_accumulated=.false.,source_funded=.false.
    ! Preserve the predictor before new-time funding changes the accepted mean.
    type(amr_tree_state_3d) :: predictor_rate,predictor_scalar_rate,predictor_pressure_work
    type(eb_integrated_flux_3d), allocatable :: predictor_faces(:),predictor_scalar_faces(:)
    type(amr_tree_state_3d) :: predictor,rate,scalar_rate,pressure_work
    type(eb_integrated_flux_3d), allocatable :: faces(:),scalar_faces(:)
    type(state_redist_register_3d), allocatable :: registers(:),scalar_registers(:)
  end type eb_godunov_hydro_stage_3d
  type, public :: eb_godunov_interval_3d
    integer :: depth=-1,iterations=0,parent_record=0,subtree_last=0
    integer :: accepted_source_stages=0
    logical :: source_prefix_only=.false.
    real(dp) :: clock_energy_correction=0,clock_pressure_work_change=0,closed_energy_correction=0
    real(dp) :: closed_pressure_work_change=0
    real(dp), allocatable :: clock_start(:),clock_interval(:)
    logical :: routed_auxiliary=.false.,replayed=.false.
    type(amr_tree_state_3d), allocatable :: clock_routed_pending(:)
    type(amr_tree_state_3d) :: clock_heat,closed_gas,closed_temperature,closed_internal_energy,closed_reaction_source
    type(amr_tree_state_3d) :: clock_reaction_source
    type(amr_tree_state_3d) :: old_predictor,new_predictor
    real(dp) :: clock_minimum_theta=1,minimum_theta=1
    ! Actual accepted off-level auxiliary impulse, after fluid/leaf filtering.
    ! Native source stages retain only the selected-level packet.
    type(amr_tree_state_3d) :: routed_auxiliary_impulse
    ! Accepted off-level old/last source rates before temporal integration.
    ! Partial common clocks cannot infer these from a selected-level packet.
    type(amr_tree_state_3d) :: routed_auxiliary_old_rate,routed_auxiliary_new_rate
    type(state_redist_register_3d), allocatable :: open_ancestor_registers(:),open_ancestor_scalar_registers(:)
    type(eb_integrated_flux_3d), allocatable :: closed_faces(:),closed_scalar_faces(:)
    real(dp) :: time=0,interval=0,face_weight=1,force(3)=0
    logical :: chemistry=.false.,transport_flags(4)=.false.,prepared_redistribution=.false.
    logical :: has_auxiliary=.false.,has_spatial_provider=.false.,spatial_diffusion_is_source=.false.
    logical :: has_valid_diffusion_provider=.false.,has_predictor_provider=.false.
    logical :: auxiliary_has_internal_energy=.false.
    logical :: has_source_funding=.false.
    type(amr_tree_state_3d) :: source_funding_stock,source_funding_impulse
    ! The wrapper clock anchor precedes target receipt payment. Restore its
    ! whole pending queue and every open ancestor family with gas/H/T.
    type(eb_held_funding_context_3d) :: clock_source_funding_context
    type(amr_patch_tree_3d) :: topology
    type(dual_energy_policy) :: policy
    type(amr_tree_state_3d) :: clock_gas,clock_temperature,clock_internal_energy,initial_reaction_source
    type(amr_tree_state_3d) :: initial_gas,initial_temperature,initial_internal_energy
    type(amr_tree_state_3d) :: final_gas,final_temperature,final_internal_energy,reaction_source
    type(amr_tree_state_3d) :: ancestor_start,ancestor_end,ancestor_temperature_start,ancestor_temperature_end
    type(amr_tree_state_3d) :: ancestor_internal_energy_start,ancestor_internal_energy_end,source_start,source_end
    real(dp), allocatable :: fractions_start(:),fractions_end(:)
    type(eb_godunov_source_stage_3d), allocatable :: sources(:)
    type(eb_godunov_hydro_stage_3d), allocatable :: hydro(:)
    type(eb_integrated_flux_3d), allocatable :: faces(:),scalar_faces(:),registered_faces(:)
    type(state_redist_register_3d), allocatable :: registers(:),scalar_registers(:)
    real(dp), allocatable :: boundary_change(:),wall_change(:),chemical_change(:),source_change(:),auxiliary_change(:)
    real(dp) :: initial_energy_correction=0,numerical_energy_change=0,pressure_work_change=0
    real(dp) :: internal_energy_projection_change=0,internal_energy_closure_change=0
    real(dp) :: initial_internal_energy_change=0
    ! Wrapper EOS impulses occur at clock entry and after all descendants.
    ! Prefixes become valid only for a completely observed accepted forest.
    real(dp) :: held_clock_start_change=0,held_clock_end_change=0
    real(dp) :: held_clock_projection_before=0,held_clock_reset_before=0
    real(dp) :: held_clock_projection_after=0,held_clock_reset_after=0
    logical :: held_clock_events_closed=.false.,held_clock_prefixes_valid=.false.
  end type eb_godunov_interval_3d
  ! Each actual prepared evaluation has a nonlinear StateRedist operator.
  ! Keep its native inputs, paired faces/registers and complete volume source;
  ! neither a final mean nor a frozen register delta can replace this family.
  type, public :: eb_mol_held_prepared_stage_3d
    integer :: source_stage=-1,depth=-1,max_order=0
    real(dp) :: interval=0,target_volume_fraction=.5_dp,register_scale=.5_dp,spatial_weight=1
    type(amr_tree_state_3d) :: state_halos,scalar_state_halos,external_rate,auxiliary_rate
    type(amr_tree_state_3d) :: raw,rate,scalar_raw,scalar_rate,low,scalar_low
    type(eb_integrated_flux_3d), allocatable :: faces(:),scalar_faces(:),low_faces(:),scalar_low_faces(:)
    type(state_redist_register_3d), allocatable :: registers(:),scalar_registers(:)
  end type eb_mol_held_prepared_stage_3d
  ! One held MOL clock. Source/face payloads describe the accepted
  ! selected level; clock payloads retain every ancestor and descendant.
  ! Prepared StateRedist retains actual stages and uses native stage replay.
  type, public :: eb_mol_held_interval_3d
    logical :: separate_momentum=.false.,consistent_velocity=.false.,paired_pressure=.false.,connected_directions=.false.
    integer :: depth=-1
    ! Declared auxiliary packets at every native MOL evaluation. A fresh
    ! local replay freezes these inputs and rebuilds force work by default.
    ! Parent source-frozen trials can also freeze external_stages explicitly;
    ! every numerical loan/face family is still rebuilt. Zero count means no cache.
    integer :: native_iterations=0
    ! Accepted ordinary native gas/UEINT spatial blend for this clock.
    real(dp) :: native_spatial_weight=1
    real(dp) :: native_force(3)=0
    type(amr_tree_state_3d), allocatable :: auxiliary_stages(:)
    ! Source-frozen native replay retains its accepted volume forcing stages;
    ! body-work and rotating momentum packets depend on the stage state.
    type(amr_tree_state_3d), allocatable :: external_stages(:)
    logical :: native_providers(2)=.false.
    integer :: native_provider_halo_width=6
    ! Use the common provider packet (rate/halos/faces/limiter), with one
    ! independent slot per spatial/direct callback at each actual stage.
    type(eb_godunov_source_stage_3d), allocatable :: provider_stages(:,:)
    real(dp) :: time=0,interval=0
    logical :: replayed=.false.,has_frame=.false.,has_auxiliary_source=.false.
    logical :: prepared_redistribution=.false.,accepted_coarse_stages=.false.
    logical :: spatial_diffusion_is_source=.false.,auxiliary_has_internal_energy=.false.,routed_auxiliary=.false.
    type(eb_mol_held_prepared_stage_3d), allocatable :: prepared_stages(:)
    type(rotating_frame_3d) :: frame
    type(dual_energy_policy) :: policy
    type(eb_reactive_tree_3d) :: initial,local_tree,clock_initial
    ! A sub-parent replay also needs the still-active ancestor trajectories;
    ! average-down of a completed descendant tail cannot reconstruct them.
    type(eb_reactive_tree_3d) :: ancestor_start,ancestor_end
    type(amr_tree_state_3d) :: ancestor_internal_energy_start,ancestor_internal_energy_end
    type(amr_tree_state_3d) :: ancestor_source_start,ancestor_source_end
    real(dp), allocatable :: fractions_start(:),fractions_end(:)
    ! Actual native clocks, including still-open ancestors. Fractions alone
    ! cannot reconstruct their periods without changing represented times.
    real(dp), allocatable :: clock_start(:),clock_interval(:)
    type(amr_tree_state_3d) :: routed_source
    type(amr_tree_state_3d) :: forcing,mass_source,initial_internal_energy,local_internal_energy
    type(amr_tree_state_3d) :: internal_energy_spatial_rate,initial_reaction_source,reaction_source
    type(amr_tree_state_3d) :: old_predictor,new_predictor
    type(amr_tree_state_3d) :: clock_internal_energy,clock_reaction_source,clock_heat
    type(state_redist_register_3d), allocatable :: open_ancestor_registers(:),open_ancestor_held_registers(:)
    real(dp) :: clock_energy_correction=0,forcing_energy_correction=0,numerical_energy_correction=0
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
    type(eb_integrated_flux_3d), allocatable :: old_flux(:),spatial_flux(:),old_internal_energy_flux(:),internal_energy_flux(:)
    real(dp), allocatable :: chemical_change(:)
    real(dp), allocatable :: wall_change(:),source_change(:),auxiliary_change(:)
  end type eb_mol_held_interval_3d
  public :: eb_auxiliary_nvar,pack_eb_held_auxiliary_state_3d,apply_eb_held_auxiliary_increment_3d
  public :: assemble_reactive_eb_held_funding_debt_3d
  public :: blend_native_mol_face_pair_3d
  public :: valid_native_mol_clock_profile_3d
  public :: reconstruct_reactive_eb_mol_held_interval_3d,replay_reactive_eb_mol_held_interval_3d
  public :: reconstruct_reactive_eb_mol_held_registration_3d
  public :: replay_reactive_eb_mol_held_group_3d
  public :: redistribute_eb_tree_state_3d,redistribute_sdc_eb_tree_state_rate_3d
  public :: redistribute_sdc_eb_tree_state_level_rate_3d
  public :: synchronize_reactive_eb_tree_register_flux_3d
  public :: build_eb_tree_extended_geometry_3d
  public :: tag_reactive_eb_tree_temperature_3d
  public :: remap_reactive_eb_tree_3d,average_down_reactive_eb_tree_field_3d
  public :: redistribute_sdc_eb_tree_flux_rate_3d,redistribute_sdc_eb_tree_level_rate_3d
  public :: reactive_eb_tree_timestep_3d
  public :: advance_reactive_eb_tree_full_3d,advance_reactive_eb_tree_sdc_3d
  public :: advance_reactive_eb_tree_sdc_subcycle_3d
  public :: advance_reactive_eb_tree_transport_3d,compute_sdc_eb_tree_transport_3d
  public :: advance_reactive_eb_tree_subcycle_transport_3d
  public :: advance_reactive_eb_tree_subcycle_hydro_3d
  public :: compute_godunov_state_redist_eb_tree_hydro_3d
  public :: advance_reactive_eb_tree_hydro_3d,compute_sdc_eb_tree_hydro_3d
  public :: compute_sdc_eb_tree_transport_halos_3d
  public :: compute_eb_tree_wall_transport_3d
  public :: compute_mol_eb_tree_hydro_3d
  public :: apply_eb_tree_flux_correction_3d
  public :: repair_reactive_eb_tree_interface_budget_3d
  public :: repair_reactive_eb_tree_held_interface_budget_3d
  public :: assemble_reactive_eb_tree_reflux_3d
  public :: synchronize_reactive_eb_tree_flux_3d
  public :: synchronize_reactive_eb_tree_scalar_flux_3d,synchronize_reactive_eb_tree_scalar_register_flux_3d
  public :: fill_reactive_eb_tree_halo_3d,fill_sdc_eb_tree_source_halo_3d
  public :: initialize_eb_tree_reaction_source_3d
  public :: react_reactive_eb_tree_3d,react_sdc_eb_tree_3d
  public :: initialize_reactive_eb_tree_3d,average_down_reactive_eb_tree_3d,integrate_reactive_eb_tree_3d
  ! Auxiliary sources receive the old state at stage 0 and provisional new
  ! states at stages 1..iterations. Rates are per fluid volume per time.
  ! Callback-owned physical state must remain private until the driver succeeds.
  ! With auxiliary_has_internal_energy enabled, the callback input state,
  ! output rate and auxiliary_change have Nv+1 rows: rows1:Nv retain the
  ! existing gas layout, and the final row is actual held UEINT/source/budget.
  ! This row is never inferred from UEDEN or interpreted as an extra species.
  abstract interface
    ! Observe the ordinary held closure before/after reflux, restriction and
    ! energy reset. Fluxes are integrated over interval, not instantaneous.
    subroutine eb_mol_sync_observer(phase,depth,time,interval,state,stored,flux,scalar_flux, &
        raw_increment,applied_increment)
      import :: dp,eb_reactive_tree_3d,amr_tree_state_3d,eb_integrated_flux_3d
      integer, intent(in) :: phase,depth
      real(dp), intent(in) :: time,interval
      type(eb_reactive_tree_3d), intent(in) :: state
      type(amr_tree_state_3d), intent(in) :: stored
      type(amr_tree_state_3d), optional, intent(in) :: raw_increment,applied_increment
      type(eb_integrated_flux_3d), intent(in) :: flux(:),scalar_flux(:)
    end subroutine
    ! Read-only attempt diagnostics, not accepted-step receipts. time is the
    ! interval start; stage zero evaluates old state, positive stages new state.
    subroutine eb_mol_stage_observer(stage,depth,time,interval,state,diagnostics,flux)
      import :: dp,eb_reactive_tree_3d,eb_mol_diagnostics_3d,eb_integrated_flux_3d
      integer, intent(in) :: stage,depth
      real(dp), intent(in) :: time,interval
      type(eb_reactive_tree_3d), intent(in) :: state
      type(eb_mol_diagnostics_3d), intent(in) :: diagnostics(:)
      type(eb_integrated_flux_3d), intent(in) :: flux(:)
    end subroutine

    ! Successful ordinary forced-reactor result, before energy recovery.
    ! Spatial sources contain hydro plus diffusion; forcing is the actual
    ! admitted rate passed to the reactor, including any reset/funding terms.
    subroutine eb_mol_reaction_observer(stage,depth,time,interval,old,updated,forcing, &
        spatial_old,spatial_new,external_old,external_new,reaction)
      import :: dp,eb_reactive_tree_3d,amr_tree_state_3d
      integer, intent(in) :: stage,depth
      real(dp), intent(in) :: time,interval
      type(eb_reactive_tree_3d), intent(in) :: old
      type(amr_tree_state_3d), intent(in) :: updated,forcing,spatial_old,spatial_new,external_old,external_new,reaction
    end subroutine

    subroutine eb_sdc_auxiliary_source(state,interval,stage,depth,rate,ok)
      import :: eb_reactive_tree_3d,amr_tree_state_3d,dp
      type(eb_reactive_tree_3d), intent(in) :: state
      real(dp), intent(in) :: interval
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      logical, intent(out) :: ok
    end subroutine
    ! Time is the level interval start relative to the outer step, for all stages.
    subroutine eb_sdc_recursive_source(state,interval,stage,depth,time,rate,ok)
      import :: eb_reactive_tree_3d,amr_tree_state_3d,dp
      type(eb_reactive_tree_3d), intent(in) :: state
      real(dp), intent(in) :: interval,time
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      logical, intent(out) :: ok
    end subroutine
    ! Callback outputs are one-based valid rates, grow6 spatial rates and
    ! integrated face fluxes. The prepared old state may contain wider donors.
    subroutine eb_sdc_spatial_diffusion_source(state,geometry,old_halos,interval,depth, &
        rate,spatial_rate,flux,limiter,ok,frame)
      import :: eb_reactive_tree_3d,eb_geometry_3d,amr_tree_state_3d,eb_integrated_flux_3d,dp,rotating_frame_3d
      type(eb_reactive_tree_3d), intent(in) :: state
      type(eb_geometry_3d), intent(in) :: geometry(:)
      type(amr_tree_state_3d), intent(in) :: old_halos
      real(dp), intent(in) :: interval
      integer, intent(in) :: depth
      type(amr_tree_state_3d), intent(inout) :: rate,spatial_rate
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
      real(dp), intent(inout) :: limiter
      logical, intent(out) :: ok
      type(rotating_frame_3d), optional, intent(in) :: frame
    end subroutine eb_sdc_spatial_diffusion_source
  end interface
  public :: eb_sdc_auxiliary_source,eb_sdc_recursive_source,eb_sdc_spatial_diffusion_source
  public :: valid_godunov_source_cache_3d,valid_godunov_clock_anchor_3d,same_godunov_predictor_3d
  public :: valid_native_mol_source_prefix_3d
  public :: rebuild_godunov_held_clock_prefixes_3d
  public :: move_godunov_source_stage_3d,extend_godunov_source_stage_3d
  public :: move_godunov_interval_3d,copy_godunov_field_3d,copy_godunov_fluxes_3d,copy_godunov_registers_3d
contains
  ! A source prefix binds observed packets to one physical clock entry.
  ! This validates packet reuse, not an accepted trajectory or wrapper anchor.
  logical function valid_native_mol_source_prefix_3d(saved,tree,geometry,depth,interval,iterations,chemistry, &
      transport_flags,force,acceleration,policy,gas,temp,held,ir,local,correction,has_auxiliary,held_auxiliary, &
      frame,funding_stock,separate_momentum,consistent_velocity,paired_pressure,connected_directions) result(good)
    logical, optional, intent(in) :: separate_momentum,consistent_velocity,paired_pressure,connected_directions
    logical :: selected_momentum,selected_velocity,selected_pressure,selected_connected
    type(eb_native_mol_source_prefix_3d), intent(in) :: saved
    type(amr_patch_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: geometry(:)
    integer, intent(in) :: depth,iterations
    real(dp), intent(in) :: interval,force(3),acceleration(3),correction
    logical, intent(in) :: chemistry,transport_flags(4),local(:),has_auxiliary,held_auxiliary
    type(dual_energy_policy), intent(in) :: policy
    type(amr_tree_state_3d), intent(in) :: gas,temp,held,ir
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d), optional, intent(in) :: funding_stock
    integer(int64), allocatable :: origin(:,:),domain(:,:),saved_origin(:,:),saved_domain(:,:)
    integer, allocatable :: ratios(:),saved_ratios(:)
    integer :: p,stage,nv,naux
    logical :: valid
    good=.false.
    selected_momentum=.false.;selected_velocity=.false.;selected_pressure=.false.;selected_connected=.false.
    if(present(separate_momentum)) selected_momentum=separate_momentum
    if(present(consistent_velocity)) selected_velocity=consistent_velocity
    if(present(paired_pressure)) selected_pressure=paired_pressure
    if(present(connected_directions)) selected_connected=connected_directions
    if(selected_velocity.and..not.selected_momentum) return
    if(selected_connected.and..not.selected_momentum) return
    if(selected_pressure.and..not.selected_velocity) return
    if(selected_pressure.and.present(frame)) return
    if(saved%separate_momentum.neqv.selected_momentum) return
    if(saved%consistent_velocity.neqv.selected_velocity) return
    if(saved%paired_pressure.neqv.selected_pressure) return
    if(saved%connected_directions.neqv.selected_connected) return
    if(.not.tree%is_valid().or..not.saved%topology%is_valid()) return
    if(size(saved%topology%nodes)/=size(tree%nodes).or.size(local)/=size(tree%nodes)) return
    if(size(geometry)/=size(local).or..not.allocated(saved%geometry)) return
    if(size(saved%geometry)/=size(local)) return
    if(saved%depth/=depth.or.saved%iterations/=iterations.or.saved%interval/=interval) return
    if(saved%accepted_source_stages<1.or.saved%accepted_source_stages>iterations+1) return
    if((saved%chemistry.neqv.chemistry).or.any(saved%transport_flags.neqv.transport_flags)) return
    if(any(saved%force/=force).or.any(saved%acceleration/=acceleration)) return
    if(.not.ieee_is_finite(correction).or.saved%initial_energy_correction/=correction) return
    if((saved%has_auxiliary_source.neqv.has_auxiliary).or. &
      (saved%auxiliary_has_internal_energy.neqv.held_auxiliary)) return
    if((saved%has_frame.neqv.present(frame)).or.(saved%has_source_funding.neqv.present(funding_stock))) return
    if(present(frame)) then
      if(saved%frame%axis/=frame%axis.or.saved%frame%angular_speed/=frame%angular_speed.or. &
        any(saved%frame%axis_point/=frame%axis_point)) return
    endif
    if(saved%policy%eta2/=policy%eta2.or.(saved%policy%allow_small_energy.neqv.policy%allow_small_energy).or. &
      (saved%policy%allow_negative_energy.neqv.policy%allow_negative_energy).or. &
      (saved%policy%update_total_from_internal.neqv.policy%update_total_from_internal)) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    if(.not.valid) return
    call amr_tree_integer_layout_3d(saved%topology,saved_origin,saved_domain,saved_ratios,valid)
    if(.not.valid) return
    if(size(ratios)/=size(saved_ratios)) return
    if(any(origin/=saved_origin).or.any(domain/=saved_domain).or.any(ratios/=saved_ratios)) return
    do p=1,size(local)
      if(tree%nodes(p)%level/=saved%topology%nodes(p)%level.or. &
        tree%nodes(p)%parent/=saved%topology%nodes(p)%parent.or. &
        any(tree%nodes(p)%shape/=saved%topology%nodes(p)%shape)) return
      if(.not.geometry(p)%is_valid().or..not.saved%geometry(p)%is_valid()) return
      if(.not.matching_eb_metrics_3d(geometry(p),saved%geometry(p))) return
      ! Source cell coordinates and the fluid mask must match exactly.
      if(any(geometry(p)%volume_fraction/=saved%geometry(p)%volume_fraction)) return
      if(any([geometry(p)%x_lower,geometry(p)%y_lower,geometry(p)%z_lower, &
        geometry(p)%dx,geometry(p)%dy,geometry(p)%dz]/= &
        [saved%geometry(p)%x_lower,saved%geometry(p)%y_lower,saved%geometry(p)%z_lower, &
        saved%geometry(p)%dx,saved%geometry(p)%dy,saved%geometry(p)%dz])) return
    enddo
    if(.not.allocated(gas%fields)) return
    if(size(gas%fields)/=size(local)) return
    ! Idle ranks have no gas component from which to infer species count.
    nv=0
    do p=1,size(local)
      if(.not.local(p)) cycle
      if(.not.allocated(gas%fields(p)%values)) return
      nv=size(gas%fields(p)%values,1);exit
    enddo
    if(.not.same_entry(saved%initial_gas,gas).or..not.same_entry(saved%initial_temperature,temp).or. &
      .not.same_entry(saved%initial_internal_energy,held).or..not.same_entry(saved%initial_reaction_source,ir)) return
    if(present(funding_stock)) then
      if(.not.same_entry(saved%source_funding_stock,funding_stock)) return
    endif
    if(.not.allocated(saved%external_stages)) return
    if(lbound(saved%external_stages,1)/=0.or.ubound(saved%external_stages,1)/=iterations) return
    if(allocated(saved%auxiliary_stages).neqv.has_auxiliary) return
    if(has_auxiliary) then
      if(lbound(saved%auxiliary_stages,1)/=0.or.ubound(saved%auxiliary_stages,1)/=iterations) return
    endif
    naux=eb_auxiliary_nvar(nv,held_auxiliary)
    do stage=0,iterations
      if(stage<saved%accepted_source_stages) then
        if(.not.valid_packet(saved%external_stages(stage),nv)) return
        if(has_auxiliary) then
          if(.not.valid_packet(saved%auxiliary_stages(stage),naux)) return
        endif
      else
        if(allocated(saved%external_stages(stage)%fields)) return
        if(has_auxiliary) then
          if(allocated(saved%auxiliary_stages(stage)%fields)) return
        endif
      endif
    enddo
    good=.true.
  contains
    logical function same_entry(a,b) result(same)
      type(amr_tree_state_3d), intent(in) :: a,b
      integer :: p
      same=.false.
      if(.not.allocated(a%fields).or..not.allocated(b%fields)) return
      if(size(a%fields)/=size(local).or.size(b%fields)/=size(local)) return
      if(lbound(a%fields,1)/=1.or.lbound(b%fields,1)/=1) return
      do p=1,size(local)
        if(.not.local(p)) cycle
        if(.not.allocated(a%fields(p)%values).or..not.allocated(b%fields(p)%values)) return
        if(any(shape(a%fields(p)%values)/=shape(b%fields(p)%values))) return
        if(any(lbound(a%fields(p)%values)/=1).or.any(lbound(b%fields(p)%values)/=1)) return
        if(.not.all(ieee_is_finite(a%fields(p)%values)).or.any(a%fields(p)%values/=b%fields(p)%values)) return
      enddo
      same=.true.
    end function same_entry
    logical function valid_packet(packet,components) result(valid)
      type(amr_tree_state_3d), intent(in) :: packet
      integer, intent(in) :: components
      integer :: p
      valid=.false.
      if(.not.allocated(packet%fields)) return
      if(size(packet%fields)/=size(local).or.lbound(packet%fields,1)/=1) return
      do p=1,size(local)
        if(.not.local(p)) cycle
        if(tree%nodes(p)%level/=depth) cycle
        if(.not.allocated(packet%fields(p)%values)) return
        if(any(shape(packet%fields(p)%values)/=[components,tree%nodes(p)%shape])) return
        if(any(lbound(packet%fields(p)%values)/=1)) return
        if(.not.all(ieee_is_finite(packet%fields(p)%values))) return
      enddo
      valid=.true.
    end function valid_packet
  end function valid_native_mol_source_prefix_3d

  subroutine move_godunov_source_stage_3d(source,destination)
    type(eb_godunov_source_stage_3d), intent(inout) :: source,destination
    destination%index=source%index
    destination%auxiliary_accepted=source%auxiliary_accepted
    destination%provider_accepted=source%provider_accepted
    destination%coarse_registers_accumulated=source%coarse_registers_accumulated
    destination%provider_kind=source%provider_kind;destination%provider_limiter=source%provider_limiter
    call move_alloc(source%provider_rate%fields,destination%provider_rate%fields)
    call move_alloc(source%provider_halos%fields,destination%provider_halos%fields)
    call move_alloc(source%gas%fields,destination%gas%fields)
    call move_alloc(source%temperature%fields,destination%temperature%fields)
    call move_alloc(source%internal_energy%fields,destination%internal_energy%fields)
    call move_alloc(source%external%fields,destination%external%fields)
    call move_alloc(source%auxiliary%fields,destination%auxiliary%fields)
    call move_alloc(source%rate%fields,destination%rate%fields)
    call move_alloc(source%scalar_rate%fields,destination%scalar_rate%fields)
    call move_alloc(source%provider_faces,destination%provider_faces)
    call move_alloc(source%faces,destination%faces)
    call move_alloc(source%scalar_faces,destination%scalar_faces)
    call move_alloc(source%registers,destination%registers)
    call move_alloc(source%scalar_registers,destination%scalar_registers)
  end subroutine move_godunov_source_stage_3d

  subroutine extend_godunov_source_stage_3d(recovered,original)
    type(eb_godunov_source_stage_3d), intent(inout) :: recovered,original
    type(eb_godunov_source_stage_3d) :: retained
    retained%auxiliary_accepted=original%auxiliary_accepted
    retained%provider_accepted=original%provider_accepted
    retained%provider_kind=original%provider_kind;retained%provider_limiter=original%provider_limiter
    if(retained%auxiliary_accepted) call move_alloc(original%auxiliary%fields,retained%auxiliary%fields)
    if(retained%provider_accepted) then
      call move_alloc(original%provider_rate%fields,retained%provider_rate%fields)
      call move_alloc(original%provider_halos%fields,retained%provider_halos%fields)
      call move_alloc(original%provider_faces,retained%provider_faces)
    endif
    call move_godunov_source_stage_3d(recovered,original)
    if(retained%auxiliary_accepted) then
      call move_alloc(retained%auxiliary%fields,original%auxiliary%fields)
      original%auxiliary_accepted=.true.
    endif
    if(retained%provider_accepted) then
      call move_alloc(retained%provider_rate%fields,original%provider_rate%fields)
      call move_alloc(retained%provider_halos%fields,original%provider_halos%fields)
      call move_alloc(retained%provider_faces,original%provider_faces)
      original%provider_accepted=.true.
      original%provider_kind=retained%provider_kind;original%provider_limiter=retained%provider_limiter
    endif
  end subroutine extend_godunov_source_stage_3d

  subroutine move_godunov_interval_3d(source,destination)
    type(eb_godunov_interval_3d), intent(inout) :: source,destination
    destination%depth=source%depth
    destination%iterations=source%iterations
    destination%accepted_source_stages=source%accepted_source_stages
    destination%source_prefix_only=source%source_prefix_only
    destination%parent_record=source%parent_record
    destination%subtree_last=source%subtree_last
    destination%clock_energy_correction=source%clock_energy_correction
    destination%clock_pressure_work_change=source%clock_pressure_work_change
    destination%closed_energy_correction=source%closed_energy_correction
    destination%closed_pressure_work_change=source%closed_pressure_work_change
    destination%routed_auxiliary=source%routed_auxiliary
    destination%replayed=source%replayed
    destination%has_source_funding=source%has_source_funding
    call move_eb_held_funding_context_3d(source%clock_source_funding_context, &
      destination%clock_source_funding_context)
    call move_alloc(source%source_funding_stock%fields,destination%source_funding_stock%fields)
    call move_alloc(source%source_funding_impulse%fields,destination%source_funding_impulse%fields)
    call move_alloc(source%clock_routed_pending,destination%clock_routed_pending)
    call move_alloc(source%clock_start,destination%clock_start)
    call move_alloc(source%clock_interval,destination%clock_interval)
    call move_alloc(source%clock_heat%fields,destination%clock_heat%fields)
    call move_alloc(source%clock_reaction_source%fields,destination%clock_reaction_source%fields)
    call move_alloc(source%old_predictor%fields,destination%old_predictor%fields)
    call move_alloc(source%new_predictor%fields,destination%new_predictor%fields)
    destination%clock_minimum_theta=source%clock_minimum_theta
    destination%minimum_theta=source%minimum_theta
    call move_alloc(source%routed_auxiliary_impulse%fields,destination%routed_auxiliary_impulse%fields)
    call move_alloc(source%routed_auxiliary_old_rate%fields,destination%routed_auxiliary_old_rate%fields)
    call move_alloc(source%routed_auxiliary_new_rate%fields,destination%routed_auxiliary_new_rate%fields)
    call move_alloc(source%closed_gas%fields,destination%closed_gas%fields)
    call move_alloc(source%closed_temperature%fields,destination%closed_temperature%fields)
    call move_alloc(source%closed_internal_energy%fields,destination%closed_internal_energy%fields)
    call move_alloc(source%closed_reaction_source%fields,destination%closed_reaction_source%fields)
    call move_alloc(source%open_ancestor_registers,destination%open_ancestor_registers)
    call move_alloc(source%open_ancestor_scalar_registers,destination%open_ancestor_scalar_registers)
    call move_alloc(source%closed_faces,destination%closed_faces)
    call move_alloc(source%closed_scalar_faces,destination%closed_scalar_faces)
    destination%time=source%time
    destination%interval=source%interval
    destination%face_weight=source%face_weight
    destination%force=source%force
    destination%chemistry=source%chemistry
    destination%transport_flags=source%transport_flags
    destination%prepared_redistribution=source%prepared_redistribution
    destination%has_auxiliary=source%has_auxiliary
    destination%has_spatial_provider=source%has_spatial_provider
    destination%spatial_diffusion_is_source=source%spatial_diffusion_is_source
    destination%has_valid_diffusion_provider=source%has_valid_diffusion_provider
    destination%has_predictor_provider=source%has_predictor_provider
    destination%auxiliary_has_internal_energy=source%auxiliary_has_internal_energy
    call move_alloc(source%topology%nodes,destination%topology%nodes)
    destination%policy=source%policy
    call move_alloc(source%clock_gas%fields,destination%clock_gas%fields)
    call move_alloc(source%clock_temperature%fields,destination%clock_temperature%fields)
    call move_alloc(source%clock_internal_energy%fields,destination%clock_internal_energy%fields)
    call move_alloc(source%initial_reaction_source%fields,destination%initial_reaction_source%fields)
    call move_alloc(source%initial_gas%fields,destination%initial_gas%fields)
    call move_alloc(source%initial_temperature%fields,destination%initial_temperature%fields)
    call move_alloc(source%initial_internal_energy%fields,destination%initial_internal_energy%fields)
    call move_alloc(source%final_gas%fields,destination%final_gas%fields)
    call move_alloc(source%final_temperature%fields,destination%final_temperature%fields)
    call move_alloc(source%final_internal_energy%fields,destination%final_internal_energy%fields)
    call move_alloc(source%reaction_source%fields,destination%reaction_source%fields)
    call move_alloc(source%ancestor_start%fields,destination%ancestor_start%fields)
    call move_alloc(source%ancestor_end%fields,destination%ancestor_end%fields)
    call move_alloc(source%ancestor_temperature_start%fields,destination%ancestor_temperature_start%fields)
    call move_alloc(source%ancestor_temperature_end%fields,destination%ancestor_temperature_end%fields)
    call move_alloc(source%ancestor_internal_energy_start%fields,destination%ancestor_internal_energy_start%fields)
    call move_alloc(source%ancestor_internal_energy_end%fields,destination%ancestor_internal_energy_end%fields)
    call move_alloc(source%source_start%fields,destination%source_start%fields)
    call move_alloc(source%source_end%fields,destination%source_end%fields)
    call move_alloc(source%fractions_start,destination%fractions_start)
    call move_alloc(source%fractions_end,destination%fractions_end)
    call move_alloc(source%sources,destination%sources)
    call move_alloc(source%hydro,destination%hydro)
    call move_alloc(source%faces,destination%faces)
    call move_alloc(source%scalar_faces,destination%scalar_faces)
    call move_alloc(source%registered_faces,destination%registered_faces)
    call move_alloc(source%registers,destination%registers)
    call move_alloc(source%scalar_registers,destination%scalar_registers)
    call move_alloc(source%boundary_change,destination%boundary_change)
    call move_alloc(source%wall_change,destination%wall_change)
    call move_alloc(source%chemical_change,destination%chemical_change)
    call move_alloc(source%source_change,destination%source_change)
    call move_alloc(source%auxiliary_change,destination%auxiliary_change)
    destination%initial_energy_correction=source%initial_energy_correction
    destination%numerical_energy_change=source%numerical_energy_change
    destination%pressure_work_change=source%pressure_work_change
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
  end subroutine move_godunov_interval_3d
  ! Side-effect-free local cache guard, shared by serial and collective MPI
  ! family preflight. Never address an absent, malformed or nonfinite packet.
  ! Compare the actual pre-corrector operator family, including sparse faces.
  ! Shape/allocation checks precede every value access; a changed boundary,
  ! ancestor halo or spatial operator must cause funding to be rebuilt.
  logical function same_godunov_predictor_3d(saved,forcing,rate,scalar_rate,work,faces,scalar_faces) result(same)
    type(eb_godunov_hydro_stage_3d), intent(in) :: saved
    type(amr_tree_state_3d), intent(in) :: forcing,rate,scalar_rate,work
    type(eb_integrated_flux_3d), allocatable, intent(in) :: faces(:),scalar_faces(:)
    same=.false.
    if(.not.field_equal(saved%predictor,forcing)) return
    if(.not.field_equal(saved%predictor_rate,rate)) return
    if(.not.field_equal(saved%predictor_scalar_rate,scalar_rate)) return
    if(.not.field_equal(saved%predictor_pressure_work,work)) return
    if(.not.faces_equal(saved%predictor_faces,faces)) return
    if(.not.faces_equal(saved%predictor_scalar_faces,scalar_faces)) return
    same=.true.
  contains
    logical function field_equal(a,b) result(equal)
      type(amr_tree_state_3d), intent(in) :: a,b
      integer :: p
      equal=.false.
      if(allocated(a%fields).neqv.allocated(b%fields)) return
      if(.not.allocated(a%fields)) then
        equal=.true.;return
      endif
      if(any(shape(a%fields)/=shape(b%fields)).or.any(lbound(a%fields)/=lbound(b%fields))) return
      do p=lbound(a%fields,1),ubound(a%fields,1)
        if(allocated(a%fields(p)%values).neqv.allocated(b%fields(p)%values)) return
        if(.not.allocated(a%fields(p)%values)) cycle
        if(any(shape(a%fields(p)%values)/=shape(b%fields(p)%values)).or. &
          any(lbound(a%fields(p)%values)/=lbound(b%fields(p)%values))) return
        if(.not.all(ieee_is_finite(a%fields(p)%values))) return
        if(any(a%fields(p)%values/=b%fields(p)%values)) return
      enddo
      equal=.true.
    end function field_equal
    logical function faces_equal(a,b) result(equal)
      type(eb_integrated_flux_3d), allocatable, intent(in) :: a(:),b(:)
      integer :: p
      equal=.false.
      if(allocated(a).neqv.allocated(b)) return
      if(.not.allocated(a)) then
        equal=.true.;return
      endif
      if(any(shape(a)/=shape(b)).or.any(lbound(a)/=lbound(b))) return
      do p=lbound(a,1),ubound(a,1)
        if(.not.face_equal(a(p)%x,b(p)%x)) return
        if(.not.face_equal(a(p)%y,b(p)%y)) return
        if(.not.face_equal(a(p)%z,b(p)%z)) return
      enddo
      equal=.true.
    end function faces_equal
    logical function face_equal(a,b) result(equal)
      real(dp), allocatable, intent(in) :: a(:,:,:,:),b(:,:,:,:)
      equal=.false.
      if(allocated(a).neqv.allocated(b)) return
      if(.not.allocated(a)) then
        equal=.true.;return
      endif
      if(any(shape(a)/=shape(b)).or.any(lbound(a)/=lbound(b))) return
      if(.not.all(ieee_is_finite(a))) return
      equal=all(a==b)
    end function face_equal
  end function same_godunov_predictor_3d

  subroutine rebuild_godunov_held_clock_prefixes_3d(records,ok)
    type(eb_godunov_interval_3d), intent(inout) :: records(:)
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
  end subroutine rebuild_godunov_held_clock_prefixes_3d

  subroutine rebuild_serial_mol_held_clock_prefixes_3d(records,ok)
    type(eb_mol_held_interval_3d), intent(inout) :: records(:)
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
  end subroutine rebuild_serial_mol_held_clock_prefixes_3d

  logical function valid_godunov_source_cache_3d(saved,tree,nv,naux,rank,owners,allow_partial) result(valid)
    type(eb_godunov_interval_3d), intent(in) :: saved
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: nv,naux
    integer, optional, intent(in) :: rank,owners(:)
    logical, optional, intent(in) :: allow_partial
    logical, allocatable :: funding_local(:)
    integer :: funding_status
    logical :: partial
    type(eb_geometry_3d) :: dimensions
    integer :: stage,p,d(3),kind,expected
    logical :: local,complete,frontier,auxiliary_known,provider_known
    valid=.false.;partial=.false.
    if(present(allow_partial)) partial=allow_partial
    if(saved%iterations<1) return
    if(.not.all(ieee_is_finite([saved%internal_energy_projection_change,saved%internal_energy_closure_change, &
      saved%initial_internal_energy_change,saved%held_clock_start_change,saved%held_clock_end_change, &
      saved%held_clock_projection_before,saved%held_clock_reset_before, &
      saved%held_clock_projection_after,saved%held_clock_reset_after]))) return
    if(saved%held_clock_prefixes_valid.and..not.saved%held_clock_events_closed) return
    if(saved%accepted_source_stages<0.or.saved%accepted_source_stages>saved%iterations+1) return
    if(.not.partial.and.(saved%source_prefix_only.or.saved%accepted_source_stages/=saved%iterations+1)) return
    if(saved%accepted_source_stages/=saved%iterations+1.and..not.saved%source_prefix_only) return
    if(present(rank).neqv.present(owners)) return
    if(.not.allocated(tree%nodes).or..not.allocated(saved%topology%nodes).or..not.allocated(saved%sources)) return
    if(size(saved%topology%nodes)/=size(tree%nodes).or.lbound(saved%topology%nodes,1)/=1) return
    if(lbound(saved%sources,1)/=0.or.ubound(saved%sources,1)/=saved%iterations) return
    if(present(owners)) then
      if(size(owners)/=size(tree%nodes)) return
    endif
    do p=1,size(tree%nodes)
      if(saved%topology%nodes(p)%parent/=tree%nodes(p)%parent.or. &
        saved%topology%nodes(p)%level/=tree%nodes(p)%level.or. &
        any(saved%topology%nodes(p)%shape/=tree%nodes(p)%shape)) return
      associate(a=>saved%topology%nodes(p)%patch,b=>tree%nodes(p)%patch)
      if(any([a%coarse_nx,a%coarse_ny,a%coarse_nz,a%coarse_i_lower,a%coarse_i_upper, &
        a%coarse_j_lower,a%coarse_j_upper,a%coarse_k_lower,a%coarse_k_upper,a%refinement_ratio]/= &
        [b%coarse_nx,b%coarse_ny,b%coarse_nz,b%coarse_i_lower,b%coarse_i_upper, &
        b%coarse_j_lower,b%coarse_j_upper,b%coarse_k_lower,b%coarse_k_upper,b%refinement_ratio])) return
      end associate
    enddo
    if(allocated(saved%clock_source_funding_context%pending).neqv. &
      allocated(saved%clock_source_funding_context%open)) return
    if(allocated(saved%clock_source_funding_context%pending)) then
      if(.not.saved%has_source_funding.or..not.allocated(saved%clock_source_funding_context%geometry)) return
      allocate(funding_local(size(tree%nodes)),source=.true.,stat=funding_status)
      if(funding_status/=0) return
      if(present(rank)) funding_local=owners==rank
      if(.not.valid_eb_held_funding_context_3d(tree,saved%clock_source_funding_context%geometry,nv, &
        funding_local,saved%clock_source_funding_context,saved%depth)) return
    else
      if(allocated(saved%clock_source_funding_context%topology%nodes).or. &
        allocated(saved%clock_source_funding_context%geometry)) return
    endif
    if(saved%has_source_funding) then
      if(.not.field_valid(saved%source_funding_stock,nv+1,0)) return
      if(.not.saved%source_prefix_only) then
        if(.not.field_valid(saved%source_funding_impulse,nv+1,0)) return
      endif
    else
      if(allocated(saved%source_funding_stock%fields).or.allocated(saved%source_funding_impulse%fields)) return
    endif
    if(saved%routed_auxiliary) then
      if(.not.saved%has_auxiliary) return
      if(saved%source_prefix_only) then
        if(saved%sources(0)%auxiliary_accepted) then
          if(.not.field_valid(saved%routed_auxiliary_old_rate,naux,0)) return
        endif
        frontier=.false.
        if(saved%accepted_source_stages<=saved%iterations) then
          frontier=saved%accepted_source_stages>0.and.saved%sources(saved%accepted_source_stages)%auxiliary_accepted
        endif
        if(saved%accepted_source_stages>1.or.frontier) then
          if(.not.field_valid(saved%routed_auxiliary_new_rate,naux,0)) return
        endif
      else
        if(.not.field_valid(saved%routed_auxiliary_impulse,naux,0)) return
      endif
    endif
    if(.not.saved%source_prefix_only) then
      if(.not.field_valid(saved%initial_gas,nv,0).or..not.field_valid(saved%initial_temperature,1,0)) return
      if(.not.field_valid(saved%initial_internal_energy,1,0)) return
      if(.not.allocated(saved%hydro)) return
      if(lbound(saved%hydro,1)/=1.or.size(saved%hydro)/=saved%iterations) return
      do stage=1,saved%iterations
        associate(packet=>saved%hydro(stage))
        if(packet%index/=stage) return
        if(.not.field_valid(packet%predictor,nv,0).or..not.field_valid(packet%rate,nv,0)) return
        if(.not.field_valid(packet%scalar_rate,1,0).or..not.field_valid(packet%pressure_work,1,0)) return
        if(.not.field_valid(packet%predictor_rate,nv,0).or. &
          .not.field_valid(packet%predictor_scalar_rate,1,0).or. &
          .not.field_valid(packet%predictor_pressure_work,1,0)) return
        if(.not.faces_valid(packet%faces,nv).or..not.faces_valid(packet%scalar_faces,1)) return
        if(.not.faces_valid(packet%predictor_faces,nv).or..not.faces_valid(packet%predictor_scalar_faces,1)) return
        end associate
      enddo
    endif
    if(saved%accepted_source_stages==0) then
      if(.not.(saved%sources(0)%auxiliary_accepted.or.saved%sources(0)%provider_accepted)) return
    endif
    do stage=0,saved%iterations
      associate(packet=>saved%sources(stage))
      complete=stage<saved%accepted_source_stages
      auxiliary_known=packet%auxiliary_accepted;provider_known=packet%provider_accepted
      frontier=partial.and.saved%source_prefix_only.and.stage==saved%accepted_source_stages.and. &
        (auxiliary_known.or.provider_known)
      if(.not.complete.and..not.frontier) then
        if(packet%index/=-1.or.auxiliary_known.or.provider_known.or.packet%provider_kind/=0) return
        if(allocated(packet%auxiliary%fields).or.allocated(packet%provider_rate%fields).or. &
          allocated(packet%provider_halos%fields).or.allocated(packet%provider_faces)) return
        cycle
      endif
      if(packet%index/=stage) return
      if(complete.and.(auxiliary_known.neqv.saved%has_auxiliary)) return
      if(auxiliary_known) then
        if(.not.saved%has_auxiliary.or..not.field_valid(packet%auxiliary,naux,0)) return
      else
        if(allocated(packet%auxiliary%fields)) return
      endif
      expected=0
      if(saved%has_spatial_provider) then
        expected=1
        if(saved%spatial_diffusion_is_source.and.stage>1.and..not.saved%prepared_redistribution) expected=0
      endif
      if(saved%has_valid_diffusion_provider) expected=2
      kind=packet%provider_kind
      if(complete.and.(provider_known.neqv.(expected>0))) return
      if(.not.provider_known) then
        if(kind/=0.or.allocated(packet%provider_rate%fields).or.allocated(packet%provider_halos%fields).or. &
          allocated(packet%provider_faces)) return
        cycle
      endif
      if(saved%has_auxiliary.and..not.auxiliary_known) return
      if(kind/=expected.or.kind==0) return
      if(.not.field_valid(packet%provider_rate,nv,0)) return
      if(kind==1) then
        if(.not.field_valid(packet%provider_halos,nv,6)) return
      endif
      if(.not.ieee_is_finite(packet%provider_limiter)) return
      if(packet%provider_limiter<0.or.packet%provider_limiter>1) return
      if(.not.allocated(packet%provider_faces)) return
      if(size(packet%provider_faces)/=size(tree%nodes).or.lbound(packet%provider_faces,1)/=1) return
      do p=1,size(tree%nodes)
        local=.true.
        if(present(owners)) local=owners(p)==rank
        if(.not.local) then
          if(allocated(packet%provider_faces(p)%x).or.allocated(packet%provider_faces(p)%y).or. &
            allocated(packet%provider_faces(p)%z)) return
          cycle
        endif
        if(saved%depth>=0.and.tree%nodes(p)%level/=saved%depth) then
          if(any(packet%provider_rate%fields(p)%values/=0)) return
          cycle
        endif
        d=tree%nodes(p)%shape;dimensions%nx=d(1);dimensions%ny=d(2);dimensions%nz=d(3)
        if(.not.valid_eb_integrated_flux_3d(packet%provider_faces(p),dimensions,nv)) return
        if(kind==1) then
          if(any(packet%provider_halos%fields(p)%values(:,7:d(1)+6,7:d(2)+6,7:d(3)+6)/= &
            packet%provider_rate%fields(p)%values)) return
        endif
      enddo
      end associate
    enddo
    valid=.true.
  contains
    logical function faces_valid(faces,components) result(good)
      type(eb_integrated_flux_3d), allocatable, intent(in) :: faces(:)
      integer, intent(in) :: components
      integer :: node,dims(3)
      logical :: needed,has_faces
      type(eb_geometry_3d) :: geometry
      good=.false.
      if(.not.allocated(faces)) return
      if(size(faces)/=size(tree%nodes).or.lbound(faces,1)/=1) return
      do node=1,size(tree%nodes)
        needed=.true.
        if(present(owners)) needed=owners(node)==rank
        has_faces=allocated(faces(node)%x).or.allocated(faces(node)%y).or.allocated(faces(node)%z)
        if(.not.needed) then
          if(has_faces) return
          cycle
        endif
        ! Inactive clocks may have no faces or actual conservative loan faces.
        if(.not.has_faces.and.saved%depth>=0.and.tree%nodes(node)%level/=saved%depth) cycle
        dims=tree%nodes(node)%shape;geometry%nx=dims(1);geometry%ny=dims(2);geometry%nz=dims(3)
        if(.not.valid_eb_integrated_flux_3d(faces(node),geometry,components)) return
      enddo
      good=.true.
    end function faces_valid
    logical function field_valid(field,components,padding) result(good)
      type(amr_tree_state_3d), intent(in) :: field
      integer, intent(in) :: components,padding
      integer :: node,dims(3)
      logical :: needed
      good=.false.
      if(.not.allocated(field%fields)) return
      if(size(field%fields)/=size(tree%nodes).or.lbound(field%fields,1)/=1) return
      do node=1,size(tree%nodes)
        needed=.true.
        if(present(owners)) needed=owners(node)==rank
        if(padding>0.and.saved%depth>=0.and.tree%nodes(node)%level/=saved%depth) needed=.false.
        if(.not.needed) then
          if(allocated(field%fields(node)%values)) return
          cycle
        endif
        if(.not.allocated(field%fields(node)%values)) return
        dims=tree%nodes(node)%shape
        if(any(shape(field%fields(node)%values)/=[components,dims+2*padding]).or. &
          any(lbound(field%fields(node)%values)/=1)) return
        if(.not.all(ieee_is_finite(field%fields(node)%values))) return
      enddo
      good=.true.
    end function field_valid
  end function valid_godunov_source_cache_3d

  logical function valid_godunov_clock_anchor_3d(saved,tree,nv,naux,levels,rank,owners) result(valid)
    type(eb_godunov_interval_3d), intent(in) :: saved
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: nv,naux,levels
    integer, optional, intent(in) :: rank,owners(:)
    integer :: node,bucket,dims(3)
    logical :: local
    valid=.false.
    if(present(rank).neqv.present(owners)) return
    if(.not.allocated(tree%nodes)) return
    if(present(owners)) then
      if(size(owners)/=size(tree%nodes)) return
    endif
    if(.not.field_valid(saved%clock_gas,nv).or..not.field_valid(saved%clock_temperature,1).or. &
      .not.field_valid(saved%clock_internal_energy,1).or..not.field_valid(saved%clock_reaction_source,nv-4)) return
    if(.not.all(ieee_is_finite([saved%clock_energy_correction,saved%clock_pressure_work_change,saved%policy%eta2]))) return
    if(.not.ieee_is_finite(saved%clock_minimum_theta)) return
    if(saved%clock_minimum_theta<0.or.saved%clock_minimum_theta>1) return
    if(.not.allocated(saved%clock_start).or..not.allocated(saved%clock_interval)) return
    if(size(saved%clock_start)/=levels.or.size(saved%clock_interval)/=levels) return
    if(lbound(saved%clock_start,1)/=1.or.lbound(saved%clock_interval,1)/=1) return
    if(.not.all(ieee_is_finite(saved%clock_start)).or..not.all(ieee_is_finite(saved%clock_interval))) return
    if(any(saved%clock_interval<=0)) return
    if(allocated(saved%clock_heat%fields)) then
      if(.not.field_valid(saved%clock_heat,1)) return
    endif
    if(saved%depth<1.or.saved%depth>=levels) return
    if(.not.field_valid(saved%ancestor_start,nv).or..not.field_valid(saved%ancestor_end,nv).or. &
      .not.field_valid(saved%ancestor_temperature_start,1).or..not.field_valid(saved%ancestor_temperature_end,1).or. &
      .not.field_valid(saved%ancestor_internal_energy_start,1).or..not.field_valid(saved%ancestor_internal_energy_end,1).or. &
      .not.field_valid(saved%source_start,nv).or..not.field_valid(saved%source_end,nv)) return
    if(.not.allocated(saved%fractions_start).or..not.allocated(saved%fractions_end)) return
    if(size(saved%fractions_start)<saved%depth.or.size(saved%fractions_end)<saved%depth) return
    if(.not.all(ieee_is_finite(saved%fractions_start)).or..not.all(ieee_is_finite(saved%fractions_end))) return
    if(any(saved%fractions_start<0).or.any(saved%fractions_start>1).or. &
      any(saved%fractions_end<0).or.any(saved%fractions_end>1)) return
    if(saved%routed_auxiliary) then
      if(.not.allocated(saved%clock_routed_pending)) return
      if(lbound(saved%clock_routed_pending,1)/=0.or.ubound(saved%clock_routed_pending,1)/=levels-1) return
      do bucket=0,levels-1
        if(.not.field_valid(saved%clock_routed_pending(bucket),naux)) return
      enddo
    endif
    if(saved%prepared_redistribution) then
      if(.not.allocated(saved%open_ancestor_registers).or..not.allocated(saved%open_ancestor_scalar_registers)) return
      if(size(saved%open_ancestor_registers)/=size(tree%nodes).or. &
        size(saved%open_ancestor_scalar_registers)/=size(tree%nodes)) return
      if(lbound(saved%open_ancestor_registers,1)/=1.or.lbound(saved%open_ancestor_scalar_registers,1)/=1) return
      do node=1,size(tree%nodes)
        local=.true.
        if(present(owners)) local=owners(node)==rank
        if(.not.local) then
          if(register_payload(saved%open_ancestor_registers(node)).or. &
            register_payload(saved%open_ancestor_scalar_registers(node))) return
          cycle
        endif
        if(.not.ieee_is_finite(saved%open_ancestor_registers(node)%scale).or. &
          .not.ieee_is_finite(saved%open_ancestor_scalar_registers(node)%scale)) return
        dims=tree%nodes(node)%shape
        ! Boundary masks depend on the supplied physical/periodic faces and
        ! are checked by the driver against its freshly initialized registers.
        if(.not.register_bounds(saved%open_ancestor_registers(node)).or. &
          .not.register_bounds(saved%open_ancestor_scalar_registers(node))) return
        if(.not.saved%open_ancestor_registers(node)%is_valid(nv,dims+12)) return
        if(.not.saved%open_ancestor_scalar_registers(node)%is_valid(1,dims+12)) return
        if(saved%open_ancestor_registers(node)%scale/=1.or.saved%open_ancestor_scalar_registers(node)%scale/=1) return
      enddo
    endif
    valid=.true.
  contains
    logical function field_valid(field,components) result(good)
      type(amr_tree_state_3d), intent(in) :: field
      integer, intent(in) :: components
      integer :: p,d(3)
      logical :: local
      good=.false.
      if(.not.allocated(field%fields)) return
      if(size(field%fields)/=size(tree%nodes).or.lbound(field%fields,1)/=1) return
      do p=1,size(tree%nodes)
        local=.true.
        if(present(owners)) local=owners(p)==rank
        if(.not.local) then
          if(allocated(field%fields(p)%values)) return
          cycle
        endif
        if(.not.allocated(field%fields(p)%values)) return
        d=tree%nodes(p)%shape
        if(any(shape(field%fields(p)%values)/=[components,d]).or.any(lbound(field%fields(p)%values)/=1)) return
        if(.not.all(ieee_is_finite(field%fields(p)%values))) return
      enddo
      good=.true.
    end function field_valid
    logical function register_bounds(reg) result(good)
      type(state_redist_register_3d), intent(in) :: reg
      good=.false.
      if(allocated(reg%coarse_delta)) then
        if(any(lbound(reg%coarse_delta)/=1)) return
      endif
      if(allocated(reg%fine_delta)) then
        if(any(lbound(reg%fine_delta)/=1)) return
      endif
      if(allocated(reg%coarse_flags)) then
        if(any(lbound(reg%coarse_flags)/=1)) return
      endif
      if(allocated(reg%fine_ghost)) then
        if(any(lbound(reg%fine_ghost)/=1)) return
      endif
      if(allocated(reg%valid_cells)) then
        if(any(lbound(reg%valid_cells)/=1)) return
      endif
      good=.true.
    end function register_bounds
    logical function register_payload(reg) result(has_payload)
      type(state_redist_register_3d), intent(in) :: reg
      has_payload=allocated(reg%coarse_delta).or.allocated(reg%fine_delta).or.allocated(reg%coarse_flags).or. &
        allocated(reg%fine_ghost).or.allocated(reg%valid_cells)
    end function register_payload
  end function valid_godunov_clock_anchor_3d

  subroutine copy_godunov_field_3d(source,destination,ok)
    type(amr_tree_state_3d), intent(in) :: source
    type(amr_tree_state_3d), intent(out) :: destination
    logical, intent(out) :: ok
    integer :: p,ios
    ok=.true.
    if(.not.allocated(source%fields)) return
    allocate(destination%fields(size(source%fields)),stat=ios)
    ok=ios==0;if(.not.ok) return
    do p=1,size(source%fields)
      if(.not.allocated(source%fields(p)%values)) cycle
      allocate(destination%fields(p)%values,source=source%fields(p)%values,stat=ios)
      ok=ios==0;if(.not.ok) return
    enddo
  end subroutine copy_godunov_field_3d
  subroutine copy_godunov_fluxes_3d(source,destination,ok)
    type(eb_integrated_flux_3d), allocatable, intent(in) :: source(:)
    type(eb_integrated_flux_3d), allocatable, intent(out) :: destination(:)
    logical, intent(out) :: ok
    integer :: p,ios
    ok=.true.
    if(.not.allocated(source)) return
    allocate(destination(size(source)),stat=ios)
    ok=ios==0;if(.not.ok) return
    do p=1,size(source)
      if(allocated(source(p)%x)) then
        allocate(destination(p)%x,source=source(p)%x,stat=ios)
        ok=ios==0;if(.not.ok) return
      endif
      if(allocated(source(p)%y)) then
        allocate(destination(p)%y,source=source(p)%y,stat=ios)
        ok=ios==0;if(.not.ok) return
      endif
      if(allocated(source(p)%z)) then
        allocate(destination(p)%z,source=source(p)%z,stat=ios)
        ok=ios==0;if(.not.ok) return
      endif
    enddo
  end subroutine copy_godunov_fluxes_3d
  subroutine copy_godunov_registers_3d(source,destination,ok)
    type(state_redist_register_3d), allocatable, intent(in) :: source(:)
    type(state_redist_register_3d), allocatable, intent(out) :: destination(:)
    logical, intent(out) :: ok
    integer :: p,ios
    ok=.true.
    if(.not.allocated(source)) return
    allocate(destination(size(source)),stat=ios)
    ok=ios==0;if(.not.ok) return
    do p=1,size(source)
      destination(p)%scale=source(p)%scale
      if(allocated(source(p)%coarse_flags)) then
        allocate(destination(p)%coarse_flags,source=source(p)%coarse_flags,stat=ios)
        ok=ios==0;if(.not.ok) return
      endif
      if(allocated(source(p)%fine_ghost)) then
        allocate(destination(p)%fine_ghost,source=source(p)%fine_ghost,stat=ios)
        ok=ios==0;if(.not.ok) return
      endif
      if(allocated(source(p)%valid_cells)) then
        allocate(destination(p)%valid_cells,source=source(p)%valid_cells,stat=ios)
        ok=ios==0;if(.not.ok) return
      endif
      if(allocated(source(p)%coarse_delta)) then
        allocate(destination(p)%coarse_delta,source=source(p)%coarse_delta,stat=ios)
        ok=ios==0;if(.not.ok) return
      endif
      if(allocated(source(p)%fine_delta)) then
        allocate(destination(p)%fine_delta,source=source(p)%fine_delta,stat=ios)
        ok=ios==0;if(.not.ok) return
      endif
    enddo
  end subroutine copy_godunov_registers_3d
  subroutine pack_godunov_tree_3d(tree,gas,temp,ok)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(out) :: gas,temp
    logical, intent(out) :: ok
    integer :: p,d(3),ios
    allocate(gas%fields(size(tree%nodes)),temp%fields(size(tree%nodes)),stat=ios)
    ok=ios==0;if(.not.ok) return
    do p=1,size(tree%nodes)
      d=tree%topology%nodes(p)%shape
      allocate(gas%fields(p)%values,source=tree%nodes(p)%state,stat=ios)
      ok=ios==0;if(.not.ok) return
      allocate(temp%fields(p)%values(1,d(1),d(2),d(3)),stat=ios)
      ok=ios==0;if(.not.ok) return
      temp%fields(p)%values(1,:,:,:)=tree%nodes(p)%temperature
    enddo
  end subroutine pack_godunov_tree_3d
  pure integer function eb_auxiliary_nvar(nvar,has_internal_energy) result(count)
    integer, intent(in) :: nvar
    logical, optional, intent(in) :: has_internal_energy
    count=nvar
    if(present(has_internal_energy)) then
      if(has_internal_energy) count=count+1
    endif
  end function eb_auxiliary_nvar

  subroutine pack_eb_held_auxiliary_state_3d(gas,stored,packed,ok)
    type(amr_tree_state_3d), intent(in) :: gas,stored
    type(amr_tree_state_3d), intent(out) :: packed
    logical, intent(out) :: ok
    integer :: p,nv,d(3),ios
    ok=.false.
    if(.not.allocated(gas%fields).or..not.allocated(stored%fields)) return
    if(size(gas%fields)/=size(stored%fields).or.lbound(gas%fields,1)/=1.or.lbound(stored%fields,1)/=1) return
    allocate(packed%fields(size(gas%fields)),stat=ios);if(ios/=0) return
    do p=1,size(gas%fields)
      if(allocated(gas%fields(p)%values).neqv.allocated(stored%fields(p)%values)) return
      if(.not.allocated(gas%fields(p)%values)) cycle
      nv=size(gas%fields(p)%values,1);d=shape(gas%fields(p)%values(1,:,:,:))
      if(any(shape(stored%fields(p)%values)/=[1,d])) return
      if(any(lbound(gas%fields(p)%values)/=1).or.any(lbound(stored%fields(p)%values)/=1)) return
      allocate(packed%fields(p)%values(nv+1,d(1),d(2),d(3)),stat=ios);if(ios/=0) return
      packed%fields(p)%values(:nv,:,:,:)=gas%fields(p)%values
      packed%fields(p)%values(nv+1:nv+1,:,:,:)=stored%fields(p)%values
    enddo
    ok=.true.
  end subroutine pack_eb_held_auxiliary_state_3d

  subroutine apply_eb_held_auxiliary_increment_3d(gas,stored,increment,ok)
    type(amr_tree_state_3d), intent(inout) :: gas,stored
    type(amr_tree_state_3d), intent(in) :: increment
    logical, intent(out) :: ok
    integer :: p,c,x,y,z,nv,d(3)
    real(dp) :: old_value,delta
    ok=.false.
    if(.not.allocated(gas%fields).or..not.allocated(stored%fields).or..not.allocated(increment%fields)) return
    if(size(gas%fields)/=size(stored%fields).or.size(gas%fields)/=size(increment%fields)) return
    if(any([lbound(gas%fields,1),lbound(stored%fields,1),lbound(increment%fields,1)]/=1)) return
    do p=1,size(gas%fields)
      if(allocated(gas%fields(p)%values).neqv.allocated(stored%fields(p)%values)) return
      if(allocated(gas%fields(p)%values).neqv.allocated(increment%fields(p)%values)) return
      if(.not.allocated(gas%fields(p)%values)) cycle
      nv=size(gas%fields(p)%values,1);d=shape(gas%fields(p)%values(1,:,:,:))
      if(any(shape(stored%fields(p)%values)/=[1,d]).or. &
        any(shape(increment%fields(p)%values)/=[nv+1,d])) return
      if(any(lbound(gas%fields(p)%values)/=1).or.any(lbound(stored%fields(p)%values)/=1).or. &
        any(lbound(increment%fields(p)%values)/=1)) return
      do z=1,d(3)
      do y=1,d(2)
      do x=1,d(1)
      do c=1,nv+1
        delta=increment%fields(p)%values(c,x,y,z)
        if(.not.ieee_is_finite(delta)) return
        if(delta==0) cycle ! Preserve untouched covered/idle payload exactly.
        if(c<=nv) then
          old_value=gas%fields(p)%values(c,x,y,z)
        else
          old_value=stored%fields(p)%values(1,x,y,z)
        endif
        if(.not.ieee_is_finite(old_value)) return
        if(delta>0) then
          if(old_value>huge(1._dp)-delta) return
        else
          if(old_value< -huge(1._dp)-delta) return
        endif
        if(c<=nv) then
          gas%fields(p)%values(c,x,y,z)=old_value+delta
        else
          stored%fields(p)%values(1,x,y,z)=old_value+delta
        endif
      enddo
      enddo
      enddo
      enddo
    enddo
    ok=.true.
  end subroutine apply_eb_held_auxiliary_increment_3d

  ! Composite flux redistribution. Covered parents do not donate or receive.
  ! Each cut leaf exchanges excess with face-connected fluid leaves, including
  ! peers and coarse/fine neighbors. Paired periodic domain faces wrap; other
  ! physical faces are closed for redistribution.
  subroutine redistribute_sdc_eb_tree_flux_rate_3d(tree,rhs,rate,ok,boundary_faces)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: rhs
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    character(len=32) :: faces(6)
    type(amr_tree_state_3d) :: candidate,active,delta
    integer(int64), allocatable :: origin(:,:),domain(:,:),scale(:)
    integer, allocatable :: ratios(:)
    integer(int64) :: cell_lo(3),cell_hi(3)
    integer :: nv,p,child,n(3),lo(3),hi(3),i,j,k,pass
    real(dp) :: kappa,volume,neighbors,total
    real(dp), allocatable :: average(:),excess(:)
    logical :: valid
    ok=.false.
    if(.not.allocated(rhs%fields)) return
    if(size(rhs%fields)<1) return
    if(.not.allocated(rhs%fields(1)%values)) return
    nv=size(rhs%fields(1)%values,1)
    if(nv<1) return
    if(.not.valid_tree(tree,nv)) return
    faces='outflow'
    if(present(boundary_faces)) faces=boundary_faces
    if(.not.valid_reactive_boundaries_3d(faces)) return
    if(.not.periodic_eb_faces_match(tree%nodes(1)%geometry,faces)) return
    if(size(rhs%fields)/=size(tree%nodes)) return
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if(.not.valid) return
    allocate(average(nv),excess(nv))
    candidate=rhs
    allocate(active%fields(size(tree%nodes)),delta%fields(size(tree%nodes)),scale(size(ratios)))
    scale=domain(1,size(ratios))/domain(1,:)
    do p=1,size(tree%nodes)
      n=tree%topology%nodes(p)%shape
      if(.not.allocated(rhs%fields(p)%values)) return
      if(any(shape(rhs%fields(p)%values)/=[nv,n])) return
      allocate(active%fields(p)%values(1,n(1),n(2),n(3)),source=1._dp)
      allocate(delta%fields(p)%values(nv,n(1),n(2),n(3)),source=0._dp)
      do child=p+1,size(tree%nodes)
        if(tree%topology%nodes(child)%parent/=p) cycle
        call bounds(tree%topology,child,lo,hi)
        active%fields(p)%values(1,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=0
      end do
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(tree%nodes(p)%geometry%volume_fraction(i,j,k)<=0) then
          candidate%fields(p)%values(:,i,j,k)=0
          active%fields(p)%values(1,i,j,k)=0
        else
          if(.not.all(ieee_is_finite(rhs%fields(p)%values(:,i,j,k)))) return
        end if
      end do
      end do
      end do
    end do
    do p=1,size(tree%nodes)
      n=tree%topology%nodes(p)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(active%fields(p)%values(1,i,j,k)==0) cycle
        kappa=tree%nodes(p)%geometry%volume_fraction(i,j,k)
        if(kappa>=1) cycle
        cell_lo=(origin(:,p)+[i-1,j-1,k-1])*scale(tree%topology%nodes(p)%level+1)
        cell_hi=cell_lo+scale(tree%topology%nodes(p)%level+1)
        volume=kappa*tree%nodes(p)%geometry%dx*tree%nodes(p)%geometry%dy*tree%nodes(p)%geometry%dz
        if(.not.ieee_is_finite(volume).or.volume<=0) return
        average=volume*rhs%fields(p)%values(:,i,j,k);neighbors=0
        do pass=1,2
          call visit_neighbors()
          if(pass==1) then
            total=volume+neighbors
            if(neighbors<=0.or..not.ieee_is_finite(total)) return
            average=average/total
            candidate%fields(p)%values(:,i,j,k)=kappa*rhs%fields(p)%values(:,i,j,k)+(1-kappa)*average
            excess=volume*(1-kappa)*(rhs%fields(p)%values(:,i,j,k)-average)/neighbors
          end if
        end do
      end do
      end do
      end do
    end do
    do p=1,size(tree%nodes)
      candidate%fields(p)%values=candidate%fields(p)%values+delta%fields(p)%values
      if(.not.all(ieee_is_finite(candidate%fields(p)%values))) return
    end do
    call move_alloc(candidate%fields,rate%fields)
    ok=.true.
  contains
    subroutine visit_neighbors()
      integer :: axis,side,q,direction,lower(3),upper(3),xyz(3),a,b,c,face_cell
      integer(int64) :: face,s,qlo(3),qhi(3),first(3),last(3)
      real(dp) :: neighbor_volume
      do axis=1,3
      do side=0,1
        if(.not.open_face(tree%nodes(p)%geometry,[i,j,k],axis,side)) cycle
        face=cell_lo(axis)
        if(side==1) face=cell_hi(axis)
        if(faces(2*axis-1)=='periodic') then
          if(side==0.and.face==0) face=domain(axis,size(ratios))
          if(side==1.and.face==domain(axis,size(ratios))) face=0
        end if
        do q=1,size(tree%nodes)
          s=scale(tree%topology%nodes(q)%level+1)
          ! A leaf can touch this face only on an integer cell boundary.
          if(modulo(face,s)/=0) cycle
          qlo=origin(:,q)*s;qhi=(origin(:,q)+tree%topology%nodes(q)%shape)*s
          if(side==0.and.(face<=qlo(axis).or.face>qhi(axis))) cycle
          if(side==1.and.(face<qlo(axis).or.face>=qhi(axis))) cycle
          first=0;last=0
          do direction=1,3
            if(direction==axis) cycle
            first(direction)=max(cell_lo(direction),qlo(direction))
            last(direction)=min(cell_hi(direction),qhi(direction))
          end do
          if(any([(last(direction)<=first(direction).and.direction/=axis,direction=1,3)])) cycle
          do direction=1,3
            if(direction==axis) cycle
            lower(direction)=int(first(direction)/s-origin(direction,q))+1
            upper(direction)=int((last(direction)-1)/s-origin(direction,q))+1
          end do
          face_cell=int(face/s-origin(axis,q))+side
          lower(axis)=face_cell;upper(axis)=face_cell
          do c=lower(3),upper(3)
          do b=lower(2),upper(2)
          do a=lower(1),upper(1)
            xyz=[a,b,c]
            if(active%fields(q)%values(1,a,b,c)==0) cycle
            if(.not.open_face(tree%nodes(q)%geometry,xyz,axis,1-side)) cycle
            neighbor_volume=tree%nodes(q)%geometry%volume_fraction(a,b,c)* &
              tree%nodes(q)%geometry%dx*tree%nodes(q)%geometry%dy*tree%nodes(q)%geometry%dz
            if(pass==1) then
              average=average+neighbor_volume*rhs%fields(q)%values(:,a,b,c)
              neighbors=neighbors+neighbor_volume
            else
              delta%fields(q)%values(:,a,b,c)=delta%fields(q)%values(:,a,b,c)+excess
            end if
          end do
          end do
          end do
        end do
      end do
      end do
    end subroutine visit_neighbors
    logical function open_face(g,cell,axis,side) result(opened)
      type(eb_geometry_3d), intent(in) :: g
      integer, intent(in) :: cell(3),axis,side
      integer :: f(3)
      f=cell;f(axis)=f(axis)-1+side
      select case(axis)
      case(1)
        opened=g%x_face_fraction(f(1),f(2),f(3))>0
      case(2)
        opened=g%y_face_fraction(f(1),f(2),f(3))>0
      case(3)
        opened=g%z_face_fraction(f(1),f(2),f(3))>0
      end select
    end function open_face
  end subroutine redistribute_sdc_eb_tree_flux_rate_3d

  ! Decompose the linear composite operator by source level. The selected
  ! level receives a rate; all other levels receive a time-integrated transfer
  ! for a caller-owned interval register. No gas state is updated here.
  subroutine redistribute_sdc_eb_tree_level_rate_3d(tree,rhs,depth,dt,rate,transfer,ok,boundary_faces)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: rhs
    integer, intent(in) :: depth
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: rate,transfer
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(amr_tree_state_3d) :: selected,result,exchange
    integer :: p,nv
    logical :: good
    ok=.false.
    if(.not.tree%topology%is_valid()) return
    if(depth<0.or.depth>=tree%topology%level_count()) return
    if(.not.ieee_is_finite(dt)) return
    if(dt<=0.or..not.allocated(rhs%fields)) return
    if(size(rhs%fields)/=size(tree%topology%nodes)) return
    if(.not.allocated(rhs%fields(1)%values)) return
    nv=size(rhs%fields(1)%values,1)
    selected=rhs
    do p=1,size(tree%topology%nodes)
      if(.not.allocated(rhs%fields(p)%values)) return
      if(any(shape(rhs%fields(p)%values)/=[nv,tree%topology%nodes(p)%shape])) return
      if(tree%topology%nodes(p)%level/=depth) selected%fields(p)%values=0
    end do
    call redistribute_sdc_eb_tree_flux_rate_3d(tree,selected,result,good,boundary_faces)
    if(.not.good) return
    exchange=result
    do p=1,size(tree%topology%nodes)
      if(tree%topology%nodes(p)%level==depth) then
        exchange%fields(p)%values=0
      else
        if(dt>1) then
          if(any(abs(exchange%fields(p)%values)>huge(1._dp)/dt)) return
        end if
        exchange%fields(p)%values=dt*exchange%fields(p)%values
        result%fields(p)%values=0
      end if
    end do
    call move_alloc(result%fields,rate%fields)
    call move_alloc(exchange%fields,transfer%fields)
    ok=.true.
  end subroutine redistribute_sdc_eb_tree_level_rate_3d

  ! S(P_l R) + D_l(U)/dt, where D_l contains conservative neighborhood
  ! corrections rooted on level l. Off-level increments are returned separately.
  subroutine redistribute_sdc_eb_tree_state_level_rate_3d(tree,rhs,depth,dt,rate,transfer,ok, &
      target_volume_fraction,boundary_faces,max_order, &
      boundary_halos,boundary_geometry)
    type(amr_tree_state_3d), optional, intent(in) :: boundary_halos
    type(eb_geometry_3d), optional, intent(in) :: boundary_geometry(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: rhs
    integer, intent(in) :: depth
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: rate,transfer
    integer, optional, intent(in) :: max_order
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: target_volume_fraction
    type(amr_tree_state_3d) :: selected,result,exchange,baseline,correction
    integer :: p,nv,i,j,k,n(3),order
    logical :: good
    ok=.false.
    if(.not.tree%topology%is_valid()) return
    if(depth<0.or.depth>=tree%topology%level_count()) return
    if(.not.ieee_is_finite(dt)) return
    if(dt<=0.or..not.allocated(rhs%fields)) return
    if(size(rhs%fields)/=size(tree%topology%nodes)) return
    if(.not.allocated(rhs%fields(1)%values)) return
    nv=size(rhs%fields(1)%values,1)
    selected=rhs
    do p=1,size(tree%topology%nodes)
      if(.not.allocated(rhs%fields(p)%values)) return
      if(any(shape(rhs%fields(p)%values)/=[nv,tree%topology%nodes(p)%shape])) return
      if(tree%topology%nodes(p)%level/=depth) selected%fields(p)%values=0
    end do
    order=0
    if(present(max_order)) order=max_order
    if(order/=0.and.order/=2.and.order/=4) return
    if(order>0) then
      if(.not.valid_tree(tree,nv)) return
      ! Limited reconstruction is nonlinear. Partition the correction of the
      ! common provisional state by neighborhood root, rather than summing
      ! separately limited state and rate reconstructions.
      baseline=rhs
      do p=1,size(tree%nodes)
        n=tree%topology%nodes(p)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          baseline%fields(p)%values(:,i,j,k)=0
          if(tree%nodes(p)%geometry%volume_fraction(i,j,k)<=0) then
            selected%fields(p)%values(:,i,j,k)=0;cycle
          endif
          if(any(.not.ieee_is_finite(tree%nodes(p)%state(:,i,j,k))).or. &
            any(.not.ieee_is_finite(rhs%fields(p)%values(:,i,j,k)))) return
          baseline%fields(p)%values(:,i,j,k)=tree%nodes(p)%state(:,i,j,k)+dt*rhs%fields(p)%values(:,i,j,k)
        enddo
        enddo
        enddo
      enddo
      call redistribute_eb_tree_state_3d(tree,baseline,correction,good,target_volume_fraction,boundary_faces, &
        correction_depth=depth,max_order=order,boundary_halos=boundary_halos,boundary_geometry=boundary_geometry)
      if(.not.good) return
      result=selected;exchange=correction
      do p=1,size(tree%nodes)
        if(tree%topology%nodes(p)%level==depth) then
          result%fields(p)%values=result%fields(p)%values+correction%fields(p)%values/dt
          exchange%fields(p)%values=0
        endif
        if(.not.all(ieee_is_finite(result%fields(p)%values)).or. &
          .not.all(ieee_is_finite(exchange%fields(p)%values))) return
      enddo
      call move_alloc(result%fields,rate%fields);call move_alloc(exchange%fields,transfer%fields)
      ok=.true.;return
    endif
    call redistribute_eb_tree_state_3d(tree,selected,result,good,target_volume_fraction,boundary_faces,max_order=max_order, &
      boundary_halos=boundary_halos,boundary_geometry=boundary_geometry)
    if(.not.good) return
    baseline=selected
    do p=1,size(tree%nodes)
      baseline%fields(p)%values=tree%nodes(p)%state
    end do
    call redistribute_eb_tree_state_3d(tree,baseline,correction,good,target_volume_fraction,boundary_faces,depth, &
      max_order=max_order, &
      boundary_halos=boundary_halos,boundary_geometry=boundary_geometry)
    if(.not.good) return
    exchange=result
    do p=1,size(tree%topology%nodes)
      if(tree%topology%nodes(p)%level==depth) then
        exchange%fields(p)%values=0
        result%fields(p)%values=result%fields(p)%values+correction%fields(p)%values/dt
      else
        if(dt>1) then
          if(any(abs(exchange%fields(p)%values)>huge(1._dp)/dt)) return
        end if
        exchange%fields(p)%values=dt*exchange%fields(p)%values+correction%fields(p)%values
        result%fields(p)%values=0
      end if
    end do
    do p=1,size(tree%nodes)
      if(.not.all(ieee_is_finite(result%fields(p)%values))) return
      if(.not.all(ieee_is_finite(exchange%fields(p)%values))) return
    end do
    call move_alloc(result%fields,rate%fields)
    call move_alloc(exchange%fields,transfer%fields)
    ok=.true.
  end subroutine redistribute_sdc_eb_tree_state_level_rate_3d

  ! Weighted StateRedist on composite fluid leaves. Large outward
  ! receivers use the single-grid weights; otherwise connected neighborhoods
  ! grow to the source cell's target physical volume with overlap accounting.
  ! correction_depth returns only the conservative displacement contributed by
  ! neighborhoods rooted at that depth, with zero hidden/solid entries.
  subroutine redistribute_eb_tree_state_3d(tree,input,result,ok,target_volume_fraction,boundary_faces, &
      correction_depth,max_order, &
      boundary_halos,boundary_geometry)
    type(amr_tree_state_3d), optional, intent(in) :: boundary_halos
    type(eb_geometry_3d), optional, intent(in) :: boundary_geometry(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: input
    type(amr_tree_state_3d), intent(inout) :: result
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: target_volume_fraction
    character(len=*), optional, intent(in) :: boundary_faces(6)
    integer, optional, intent(in) :: max_order
    integer, optional, intent(in) :: correction_depth
    type :: neighborhood
      integer, allocatable :: member(:)
    end type
    type(neighborhood), allocatable :: groups(:)
    type(amr_tree_state_3d) :: ids,candidate
    type(eb_geometry_3d) :: checked_geometry
    integer(int64), allocatable :: origin(:,:),domain(:,:),scale(:)
    integer, allocatable :: ratios(:),cells(:,:),receiver(:),memberships(:),queue(:),seen(:),neighbors(:),directions(:)
    real(dp), allocatable :: volume(:),target(:),beta(:),self(:),value(:,:),average(:,:),output(:,:)
    real(dp), allocatable :: center(:,:),lattice(:,:),centroid(:,:),width(:,:),fraction_cell(:), &
      slopes(:,:,:),stencil_centroid(:,:),stencil_center(:,:),stencil_width(:,:),positions(:,:),lower_bound(:),upper_bound(:)
    real(dp), allocatable :: stencil_value(:,:),stencil_fraction(:),stencil_fit(:,:)
    real(dp), allocatable :: ghost_value(:,:)
    real(dp) :: ghost_point(3,125),ghost_fit(3,125),ghost_center(3,125),ghost_fraction(125),root_centroid(3)
    logical :: ghost_near(125),extdir(6)
    integer :: ghost_count,global_cell(3),global_dims(3),extdir_bounds(2,3),pool_size
    integer :: selected_order
    integer(int64) :: maximum,cell_lo(3),cell_hi(3)
    character(len=32) :: faces(6)
    integer :: nv,p,n(3),i,j,k,child,lo(3),hi(3),nc,root,other,m,count,head,tail,current,allocation_status
    real(dp) :: fraction,total,weight,normal(3),score,best
    logical, allocatable :: stencil_connected(:)
    logical :: valid,fallback
    ok=.false.;fraction=.5_dp;faces='outflow'
    selected_order=0
    if(present(max_order)) selected_order=max_order
    if(selected_order/=0.and.selected_order/=2.and.selected_order/=4) return
    if(present(target_volume_fraction)) fraction=target_volume_fraction
    if(present(boundary_faces)) faces=boundary_faces
    if(.not.ieee_is_finite(fraction)) return
    if(fraction<=0.or.fraction>1.or..not.valid_reactive_boundaries_3d(faces)) return
    if(.not.allocated(input%fields)) return
    if(size(input%fields)<1) return
    if(.not.allocated(input%fields(1)%values)) return
    nv=size(input%fields(1)%values,1)
    if(nv<1.or..not.valid_tree(tree,nv)) return
    if(size(input%fields)/=size(tree%nodes)) return
    if(present(correction_depth)) then
      if(correction_depth<0.or.correction_depth>=tree%topology%level_count()) return
    end if
    if(.not.periodic_eb_faces_match(tree%nodes(1)%geometry,faces)) return
    if(present(boundary_halos).neqv.present(boundary_geometry)) return
    if(present(boundary_halos)) then
      if(.not.allocated(boundary_halos%fields)) return
      if(size(boundary_halos%fields)/=size(tree%nodes).or.size(boundary_geometry)/=size(tree%nodes)) return
      do p=1,size(tree%nodes)
        if(.not.allocated(boundary_halos%fields(p)%values)) return
        n=tree%topology%nodes(p)%shape
        if(any(shape(boundary_halos%fields(p)%values)/=[nv,n+4])) return
        if(any([boundary_geometry(p)%nx,boundary_geometry(p)%ny,boundary_geometry(p)%nz]/=n+4)) return
        if(any(.not.ieee_is_finite(boundary_halos%fields(p)%values))) return
        call extract_eb_geometry_3d(boundary_geometry(p),[3,3,3],n+2,checked_geometry,valid)
        if(.not.valid) return
        if(.not.matching_eb_metrics_3d(tree%nodes(p)%geometry,checked_geometry)) return
      enddo
    endif
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if(.not.valid) return
    maximum=0
    do p=1,size(tree%nodes)
      if(.not.allocated(input%fields(p)%values)) return
      n=tree%topology%nodes(p)%shape
      if(any(shape(input%fields(p)%values)/=[nv,n])) return
      if(int(n(1),int64)*n(2)>int(huge(1),int64)/n(3)) return
      maximum=maximum+int(n(1),int64)*n(2)*n(3)
      if(maximum>huge(1)) return
    end do
    allocate(ids%fields(size(tree%nodes)),scale(size(ratios)),cells(4,int(maximum)),stat=allocation_status)
    if(allocation_status/=0) return
    scale=domain(1,size(ratios))/domain(1,:)
    candidate=input;nc=0
    if(present(correction_depth)) then
      do p=1,size(tree%nodes)
        candidate%fields(p)%values=0
      end do
    end if
    do p=1,size(tree%nodes)
      n=tree%topology%nodes(p)%shape
      allocate(ids%fields(p)%values(1,n(1),n(2),n(3)),source=0._dp,stat=allocation_status)
      if(allocation_status/=0) return
      do child=p+1,size(tree%nodes)
        if(tree%topology%nodes(child)%parent/=p) cycle
        call bounds(tree%topology,child,lo,hi)
        ids%fields(p)%values(1,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=-1
      end do
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(tree%nodes(p)%geometry%volume_fraction(i,j,k)<=0) then
          candidate%fields(p)%values(:,i,j,k)=0
          ids%fields(p)%values(1,i,j,k)=0
          cycle
        end if
        if(.not.all(ieee_is_finite(input%fields(p)%values(:,i,j,k)))) return
        if(ids%fields(p)%values(1,i,j,k)<0) cycle
        nc=nc+1;cells(:,nc)=[p,i,j,k];ids%fields(p)%values(1,i,j,k)=real(nc,dp)
      end do
      end do
      end do
    end do
    if(nc>huge(1)/6) return
    allocate(groups(nc),receiver(nc),memberships(nc),queue(nc),seen(nc),neighbors(6*nc),directions(6*nc), &
      volume(nc),target(nc),beta(nc),self(nc),value(nv,nc),average(nv,nc),output(nv,nc),stat=allocation_status)
    if(allocation_status/=0) return
    do root=1,nc
      p=cells(1,root);i=cells(2,root);j=cells(3,root);k=cells(4,root)
      associate(g=>tree%nodes(p)%geometry)
        target(root)=fraction*g%dx*g%dy*g%dz
        volume(root)=g%volume_fraction(i,j,k)*g%dx*g%dy*g%dz
      end associate
      if(min(volume(root),target(root))<=0) return
      value(:,root)=input%fields(p)%values(:,i,j,k)
    end do
    receiver=0;memberships=1;self=1;beta=0;fallback=.false.
    do root=1,nc
      if(volume(root)>=target(root)) cycle
      current=root
      call adjacent(count)
      p=cells(1,root);i=cells(2,root);j=cells(3,root);k=cells(4,root)
      associate(g=>tree%nodes(p)%geometry)
        normal=[g%boundary_normal_x(i,j,k),g%boundary_normal_y(i,j,k),g%boundary_normal_z(i,j,k)]
      end associate
      best=-1
      do m=1,count
        other=neighbors(m)
        if(volume(other)<max(target(root),target(other))) cycle
        score=normal(abs(directions(m)))*real(sign(1,directions(m)),dp)
        if(score<=512*epsilon(1._dp).or.score<=best) cycle
        best=score;receiver(root)=other
      end do
      if(receiver(root)==0) then
        fallback=.true.;exit
      end if
      other=receiver(root)
      beta(root)=(target(root)-volume(root))/volume(other)
      memberships(other)=memberships(other)+1
    end do
    if(.not.fallback) then
      do root=1,nc
        other=receiver(root)
        if(other==0) cycle
        self(other)=self(other)-beta(root)/real(memberships(other),dp)
      end do
      do root=1,nc
        total=self(root)*volume(root);average(:,root)=total*value(:,root)
        other=receiver(root)
        if(other>0) then
          weight=beta(root)*volume(other)/real(memberships(other),dp)
          total=total+weight;average(:,root)=average(:,root)+weight*value(:,other)
        end if
        if(total<=0.or.self(root)<0) return
        average(:,root)=average(:,root)/total
      end do
      output=0
      do root=1,nc
        if(present(correction_depth)) then
          if(tree%topology%nodes(cells(1,root))%level/=correction_depth) cycle
          output(:,root)=output(:,root)+self(root)*real(memberships(root),dp)*(average(:,root)-value(:,root))
          other=receiver(root)
          if(other>0) output(:,other)=output(:,other)+beta(root)*(average(:,root)-value(:,other))
        else
          output(:,root)=output(:,root)+self(root)*real(memberships(root),dp)*average(:,root)
          other=receiver(root)
          if(other>0) output(:,other)=output(:,other)+beta(root)*average(:,root)
        end if
      end do
      do root=1,nc
        output(:,root)=output(:,root)/real(memberships(root),dp)
      end do
    else
      memberships=1;seen=0
      do root=1,nc
        if(volume(root)>=target(root)) cycle
        head=1;tail=1;queue(1)=root;seen(root)=root;total=volume(root)
        do while(total<target(root).and.head<=tail)
          current=queue(head)
          call adjacent(count)
          do m=1,count
            if(total>=target(root)) exit
            other=neighbors(m)
            if(seen(other)==root) cycle
            tail=tail+1;queue(tail)=other;seen(other)=root;total=total+volume(other)
          end do
          head=head+1
        end do
        if(total<target(root)) return
        allocate(groups(root)%member(tail),stat=allocation_status)
        if(allocation_status/=0) return
        groups(root)%member=queue(:tail)
        do m=2,tail
          other=queue(m);memberships(other)=memberships(other)+1
        end do
      end do
      output=0
      do root=1,nc
        count=1
        if(allocated(groups(root)%member)) count=size(groups(root)%member)
        total=0;average(:,root)=0
        do m=1,count
          other=root
          if(allocated(groups(root)%member)) other=groups(root)%member(m)
          weight=volume(other)/real(memberships(other),dp)
          total=total+weight;average(:,root)=average(:,root)+weight*value(:,other)
        end do
        if(total<=0) return
        average(:,root)=average(:,root)/total
        if(present(correction_depth)) then
          if(tree%topology%nodes(cells(1,root))%level/=correction_depth) cycle
        end if
        do m=1,count
          other=root
          if(allocated(groups(root)%member)) other=groups(root)%member(m)
          if(present(correction_depth)) then
            output(:,other)=output(:,other)+(average(:,root)-value(:,other))/real(memberships(other),dp)
          else
            output(:,other)=output(:,other)+average(:,root)/real(memberships(other),dp)
          end if
        end do
      end do
    end if
    if(selected_order>0) then
      call reconstruct_graph_state(valid)
      if(.not.valid) return
    endif
    if(.not.all(ieee_is_finite(output))) return
    do root=1,nc
      p=cells(1,root);i=cells(2,root);j=cells(3,root);k=cells(4,root)
      candidate%fields(p)%values(:,i,j,k)=output(:,root)
    end do
    call move_alloc(candidate%fields,result%fields);ok=.true.
  contains
    ! Periodic image selection is shared by the neighborhood moment and its
    ! reconstructed recipients, so their weighted displacement sums to zero.
    function displacement(other,root) result(delta)
      integer, intent(in) :: other,root
      real(dp) :: delta(3),length
      integer :: a
      delta=center(:,other)-center(:,root)
      do a=1,3
        if(faces(2*a-1)/='periodic') cycle
        length=real(domain(a,size(ratios)),dp)
        delta(a)=modulo(delta(a)+.5_dp*length,length)-.5_dp*length
      enddo
    end function displacement
    subroutine reconstruct_graph_state(success)
      logical, intent(out) :: success
      integer :: r,t,nmember,h,last,next,c,member_count
      real(dp) :: w,sum_weight,shift(3),offset(3),raw_offset(3)
      logical :: reconstructed
      success=.false.
      if(nc>huge(1)-125) return
      pool_size=nc+125
      allocate(ghost_value(nv,125))
      allocate(stencil_value(nv,pool_size),stencil_fraction(pool_size),stencil_fit(3,pool_size))
      allocate(center(3,nc),lattice(3,nc),centroid(3,nc),width(3,nc),fraction_cell(nc),slopes(3,nv,nc), &
        stencil_centroid(3,pool_size),stencil_center(3,pool_size),stencil_width(3,pool_size), &
        positions(3,nc),stencil_connected(pool_size), &
        lower_bound(nv),upper_bound(nv),stat=allocation_status)
      if(allocation_status/=0) return
      do r=1,nc
        p=cells(1,r);i=cells(2,r);j=cells(3,r);k=cells(4,r)
        w=real(scale(tree%topology%nodes(p)%level+1),dp)
        width(:,r)=w;lattice(:,r)=(real(origin(:,p),dp)+real([i,j,k],dp)-.5_dp)*w
        associate(g=>tree%nodes(p)%geometry)
          center(:,r)=lattice(:,r)+w*[g%cell_centroid_x(i,j,k),g%cell_centroid_y(i,j,k),g%cell_centroid_z(i,j,k)]
          fraction_cell(r)=g%volume_fraction(i,j,k)
        end associate
      enddo
      do r=1,nc
        shift=0;sum_weight=0
        if(fallback) then
          nmember=1
          if(allocated(groups(r)%member)) nmember=size(groups(r)%member)
          do t=1,nmember
            next=r
            if(allocated(groups(r)%member)) next=groups(r)%member(t)
            w=volume(next)/real(memberships(next),dp)
            shift=shift+w*displacement(next,r);sum_weight=sum_weight+w
          enddo
        else
          sum_weight=self(r)*volume(r);next=receiver(r)
          if(next>0) then
            w=beta(r)*volume(next)/real(memberships(next),dp)
            shift=w*displacement(next,r);sum_weight=sum_weight+w
          endif
        endif
        if(sum_weight<=0) return
        centroid(:,r)=center(:,r)+shift/sum_weight
      enddo
      do c=1,nv
        lower_bound(c)=minval(value(c,:));upper_bound(c)=maxval(value(c,:))
      enddo
      seen=0;slopes=0
      do r=1,nc
        ! Only faces joining real fluid leaves enter this stencil. Covered
        ! parent values never become neighbors, including at refinement seams.
        h=1;last=1;queue(1)=r;seen(r)=r
        do while(h<=last)
          current=queue(h);call adjacent(member_count)
          do t=1,member_count
            next=neighbors(t)
            if(seen(next)==r) cycle
            raw_offset=displacement(next,r)-(center(:,next)-lattice(:,next))+(center(:,r)-lattice(:,r))
            if(any(abs(raw_offset/width(:,r))>2+512*epsilon(1._dp))) cycle
            last=last+1;queue(last)=next;seen(next)=r
          enddo
          h=h+1
        enddo
        stencil_connected(:last)=.true.
        ! Include noncovered leaf values for the grown AMReX stencil even
        ! when an aperture path does not connect them to the root.
        do next=1,nc
          if(seen(next)==r) cycle
          raw_offset=displacement(next,r)-(center(:,next)-lattice(:,next))+(center(:,r)-lattice(:,r))
          if(any(abs(raw_offset/width(:,r))>2+512*epsilon(1._dp))) cycle
          last=last+1;queue(last)=next;stencil_connected(last)=.false.
        enddo
        do t=1,last
          next=queue(t);offset=displacement(next,r)
          stencil_centroid(:,t)=(offset+centroid(:,next)-center(:,next))/width(:,r)
          stencil_center(:,t)=(offset-(center(:,next)-lattice(:,next))+(center(:,r)-lattice(:,r)))/width(:,r)
          stencil_width(:,t)=width(:,next)/width(:,r)
          stencil_value(:,t)=average(:,next);stencil_fraction(t)=fraction_cell(next)
        enddo
        stencil_fit(:,:last)=stencil_centroid(:,:last)
        extdir=.false.;extdir_bounds=0
        if(present(boundary_halos)) then
          p=cells(1,r);global_cell=int(origin(:,p))+cells(2:4,r)
          global_dims=int(domain(:,tree%topology%nodes(p)%level+1))
          extdir=faces=='inflow';extdir_bounds(1,:)=1-global_cell;extdir_bounds(2,:)=global_dims-global_cell
          call state_redistribution_boundary_nodes_3d(boundary_geometry(p),boundary_halos%fields(p)%values, &
            cells(2:4,r)+2,global_cell,global_dims,faces,ghost_value,ghost_point,ghost_fit,ghost_center, &
            ghost_fraction,ghost_near,ghost_count,reconstructed)
          if(.not.reconstructed) return
          root_centroid=(center(:,r)-lattice(:,r))/width(:,r)
          do t=1,ghost_count
            last=last+1;stencil_value(:,last)=ghost_value(:,t)
            stencil_fraction(last)=ghost_fraction(t);stencil_connected(last)=ghost_near(t)
            stencil_centroid(:,last)=ghost_point(:,t)-root_centroid
            stencil_fit(:,last)=ghost_fit(:,t)-root_centroid
            stencil_center(:,last)=ghost_center(:,t);stencil_width(:,last)=1
          enddo
        endif
        call state_redistribution_graph_slopes_3d(stencil_value(:,:last),stencil_centroid(:,:last), &
          stencil_center(:,:last),stencil_width(:,:last),stencil_fraction(:last),selected_order,slopes(:,:,r),reconstructed, &
          near_connected=stencil_connected(:last),extdir_faces=extdir,extdir_bounds=extdir_bounds, &
          fit_centroid=stencil_fit(:,:last))
        if(.not.reconstructed) return
        if(fallback) then
          nmember=1
          if(allocated(groups(r)%member)) nmember=size(groups(r)%member)
          do t=1,nmember
            next=r
            if(allocated(groups(r)%member)) next=groups(r)%member(t)
            positions(:,t)=displacement(next,r)/width(:,r)
          enddo
        else
          nmember=1;positions(:,1)=0;next=receiver(r)
          if(next>0) then
            nmember=2;positions(:,2)=displacement(next,r)/width(:,r)
          endif
        endif
        shift=(centroid(:,r)-center(:,r))/width(:,r)
        call limit_state_recipient_slopes(average(:,r),shift,slopes(:,:,r),positions(:,:nmember),lower_bound,upper_bound)
        slopes(:,:,r)=slopes(:,:,r)/spread(width(:,r),2,nv)
      enddo
      output=0
      do r=1,nc
        if(present(correction_depth)) then
          if(tree%topology%nodes(cells(1,r))%level/=correction_depth) cycle
        endif
        if(fallback) then
          nmember=1
          if(allocated(groups(r)%member)) nmember=size(groups(r)%member)
        else
          nmember=1
          if(receiver(r)>0) nmember=2
        endif
        do t=1,nmember
          next=r;w=1._dp/real(memberships(r),dp)
          if(fallback) then
            if(allocated(groups(r)%member)) next=groups(r)%member(t)
            w=1._dp/real(memberships(next),dp)
          else
            w=self(r)
            if(t==2) then
              next=receiver(r);w=beta(r)/real(memberships(next),dp)
            endif
          endif
          offset=displacement(next,r)-(centroid(:,r)-center(:,r))
          do c=1,nv
            sum_weight=average(c,r)+dot_product(offset,slopes(:,c,r))
            if(present(correction_depth)) sum_weight=sum_weight-value(c,next)
            output(c,next)=output(c,next)+w*sum_weight
          enddo
        enddo
      enddo
      success=all(ieee_is_finite(output))
    end subroutine reconstruct_graph_state
    subroutine adjacent(count)
      integer, intent(out) :: count
      integer :: axis,side,q,direction,lower(3),upper(3),xyz(3),a,b,c,face_cell
      integer(int64) :: face,s,qlo(3),qhi(3),first(3),last(3)
      integer :: neighbor
      p=cells(1,current);i=cells(2,current);j=cells(3,current);k=cells(4,current)
      cell_lo=(origin(:,p)+[i-1,j-1,k-1])*scale(tree%topology%nodes(p)%level+1)
      cell_hi=cell_lo+scale(tree%topology%nodes(p)%level+1)
      count=0
      do axis=1,3
      do side=0,1
        if(.not.open_face(tree%nodes(p)%geometry,[i,j,k],axis,side)) cycle
        face=cell_lo(axis)
        if(side==1) face=cell_hi(axis)
        if(faces(2*axis-1)=='periodic') then
          if(side==0.and.face==0) face=domain(axis,size(ratios))
          if(side==1.and.face==domain(axis,size(ratios))) face=0
        end if
        do q=1,size(tree%nodes)
          s=scale(tree%topology%nodes(q)%level+1)
          ! A leaf can touch this face only on an integer cell boundary.
          if(modulo(face,s)/=0) cycle
          qlo=origin(:,q)*s;qhi=(origin(:,q)+tree%topology%nodes(q)%shape)*s
          if(side==0.and.(face<=qlo(axis).or.face>qhi(axis))) cycle
          if(side==1.and.(face<qlo(axis).or.face>=qhi(axis))) cycle
          first=0;last=0
          do direction=1,3
            if(direction==axis) cycle
            first(direction)=max(cell_lo(direction),qlo(direction))
            last(direction)=min(cell_hi(direction),qhi(direction))
          end do
          if(any([(last(direction)<=first(direction).and.direction/=axis,direction=1,3)])) cycle
          do direction=1,3
            if(direction==axis) cycle
            lower(direction)=int(first(direction)/s-origin(direction,q))+1
            upper(direction)=int((last(direction)-1)/s-origin(direction,q))+1
          end do
          face_cell=int(face/s-origin(axis,q))+side
          lower(axis)=face_cell;upper(axis)=face_cell
          do c=lower(3),upper(3)
          do b=lower(2),upper(2)
          do a=lower(1),upper(1)
            xyz=[a,b,c]
            if(ids%fields(q)%values(1,a,b,c)<=0) cycle
            if(.not.open_face(tree%nodes(q)%geometry,xyz,axis,1-side)) cycle
            neighbor=int(ids%fields(q)%values(1,a,b,c))
            count=count+1;neighbors(count)=neighbor;directions(count)=axis*(2*side-1)
          end do
          end do
          end do
        end do
      end do
      end do
    end subroutine adjacent
    logical function open_face(g,cell,axis,side) result(opened)
      type(eb_geometry_3d), intent(in) :: g
      integer, intent(in) :: cell(3),axis,side
      integer :: f(3)
      f=cell;f(axis)=f(axis)-1+side
      select case(axis)
      case(1)
        opened=g%x_face_fraction(f(1),f(2),f(3))>0
      case(2)
        opened=g%y_face_fraction(f(1),f(2),f(3))>0
      case(3)
        opened=g%z_face_fraction(f(1),f(2),f(3))>0
      end select
    end function open_face
  end subroutine redistribute_eb_tree_state_3d

  ! Match the single-grid StateRedist source convention: S(U+dt*R)-U,
  ! divided by dt. Covered cells are excluded before arithmetic.
  subroutine redistribute_sdc_eb_tree_state_rate_3d(tree,rhs,dt,rate,ok,target_volume_fraction,boundary_faces,max_order, &
      boundary_halos,boundary_geometry)
    type(amr_tree_state_3d), optional, intent(in) :: boundary_halos
    type(eb_geometry_3d), optional, intent(in) :: boundary_geometry(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: rhs
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: rate
    integer, optional, intent(in) :: max_order
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: target_volume_fraction
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(amr_tree_state_3d) :: provisional,candidate
    integer :: nv,p,n(3),i,j,k
    logical :: valid
    ok=.false.
    if(.not.ieee_is_finite(dt)) return
    if(dt<=0.or..not.allocated(rhs%fields)) return
    if(size(rhs%fields)<1) return
    if(.not.allocated(rhs%fields(1)%values)) return
    nv=size(rhs%fields(1)%values,1)
    if(nv<1.or..not.valid_tree(tree,nv)) return
    if(size(rhs%fields)/=size(tree%nodes)) return
    provisional=rhs
    do p=1,size(tree%nodes)
      n=tree%topology%nodes(p)%shape
      if(.not.allocated(rhs%fields(p)%values)) return
      if(any(shape(rhs%fields(p)%values)/=[nv,n])) return
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        provisional%fields(p)%values(:,i,j,k)=0
        if(tree%nodes(p)%geometry%volume_fraction(i,j,k)<=0) cycle
        if(.not.all(ieee_is_finite(rhs%fields(p)%values(:,i,j,k)))) return
        if(.not.all(ieee_is_finite(tree%nodes(p)%state(:,i,j,k)))) return
        provisional%fields(p)%values(:,i,j,k)=tree%nodes(p)%state(:,i,j,k)+dt*rhs%fields(p)%values(:,i,j,k)
      end do
      end do
      end do
    end do
    call redistribute_eb_tree_state_3d(tree,provisional,candidate,valid,target_volume_fraction,boundary_faces, &
      max_order=max_order, &
      boundary_halos=boundary_halos,boundary_geometry=boundary_geometry)
    if(.not.valid) return
    do p=1,size(tree%nodes)
      n=tree%topology%nodes(p)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(tree%nodes(p)%geometry%volume_fraction(i,j,k)<=0) cycle
        candidate%fields(p)%values(:,i,j,k)=(candidate%fields(p)%values(:,i,j,k)-tree%nodes(p)%state(:,i,j,k))/dt
      end do
      end do
      end do
      if(.not.all(ieee_is_finite(candidate%fields(p)%values))) return
    end do
    call move_alloc(candidate%fields,rate%fields);ok=.true.
  end subroutine redistribute_sdc_eb_tree_state_rate_3d

  subroutine bounds(tree,node,lo,hi)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: node
    integer, intent(out) :: lo(3),hi(3)
    lo=[tree%nodes(node)%patch%coarse_i_lower,tree%nodes(node)%patch%coarse_j_lower, &
      tree%nodes(node)%patch%coarse_k_lower]
    hi=[tree%nodes(node)%patch%coarse_i_upper,tree%nodes(node)%patch%coarse_j_upper, &
      tree%nodes(node)%patch%coarse_k_upper]
  end subroutine bounds

  ! Build target-resolution margins from a finest domain padded by width times
  ! the total refinement ratio. Geometry remains replicated for MPI callers.
  subroutine build_eb_tree_extended_geometry_3d(topology,physical,finest,patches,ok,halo_width)
    type(amr_patch_tree_3d), intent(in) :: topology
    type(eb_geometry_3d), intent(in) :: physical(:),finest
    type(eb_geometry_3d), allocatable, intent(inout) :: patches(:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: halo_width
    type(eb_geometry_3d), allocatable :: levels(:),candidate(:)
    type(eb_geometry_3d) :: center
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer(int64) :: padding
    integer, allocatable :: ratios(:),margin(:)
    integer :: width,count,level,node,n(3),lo(3),hi(3)
    logical :: valid
    ok=.false.;width=6
    if(present(halo_width)) width=halo_width
    if(width<1.or.width>huge(1)/4) return
    call amr_tree_integer_layout_3d(topology,origin,domain,ratios,valid)
    if(.not.valid.or..not.finest%is_valid()) return
    if(size(physical)/=size(topology%nodes)) return
    count=size(ratios)
    padding=domain(1,count)/domain(1,1)
    if(padding>int(huge(1)/2,int64)/width) return
    padding=padding*width
    if(any(domain(:,count)>int(huge(1),int64)-2*padding)) return
    if(any(int([finest%nx,finest%ny,finest%nz],int64)/=domain(:,count)+2*padding)) return
    allocate(levels(count),margin(count),candidate(size(physical)))
    levels(count)=finest;margin(count)=int(padding)
    do level=count,2,-1
      if(mod(margin(level),ratios(level))/=0) return
      margin(level-1)=margin(level)/ratios(level)
      call coarsen_eb_geometry_3d(levels(level),ratios(level),levels(level-1),valid)
      if(.not.valid) return
    end do
    do node=1,size(physical)
      if(.not.physical(node)%is_valid()) return
      level=topology%nodes(node)%level+1;n=topology%nodes(node)%shape
      if(any([physical(node)%nx,physical(node)%ny,physical(node)%nz]/=n)) return
      lo=int(origin(:,node))+1+margin(level)-width;hi=lo+n+2*width-1
      call extract_eb_geometry_3d(levels(level),lo,hi,candidate(node),valid)
      if(.not.valid) return
      call extract_eb_geometry_3d(candidate(node),[width+1,width+1,width+1],n+width,center,valid)
      if(.not.valid) return
      if(.not.matching_eb_metrics_3d(center,physical(node))) return
      ! Coordinates used during clipping can differ by roundoff after padding.
      ! Keep the accepted physical patch metrics bit-for-bit, including its faces.
      associate(g=>candidate(node),p=>physical(node))
      g%volume_fraction(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%volume_fraction
      g%cell_centroid_x(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%cell_centroid_x
      g%cell_centroid_y(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%cell_centroid_y
      g%cell_centroid_z(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%cell_centroid_z
      g%boundary_area(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%boundary_area
      g%boundary_centroid_x(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%boundary_centroid_x
      g%boundary_centroid_y(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%boundary_centroid_y
      g%boundary_centroid_z(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%boundary_centroid_z
      g%boundary_normal_x(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%boundary_normal_x
      g%boundary_normal_y(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%boundary_normal_y
      g%boundary_normal_z(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%boundary_normal_z
      g%boundary_normal_integral_x(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%boundary_normal_integral_x
      g%boundary_normal_integral_y(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%boundary_normal_integral_y
      g%boundary_normal_integral_z(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%boundary_normal_integral_z
      g%cell_type(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%cell_type
      g%x_face_fraction(width:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%x_face_fraction
      g%x_face_centroid_y(width:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%x_face_centroid_y
      g%x_face_centroid_z(width:width+n(1),width+1:width+n(2),width+1:width+n(3))=p%x_face_centroid_z
      g%y_face_fraction(width+1:width+n(1),width:width+n(2),width+1:width+n(3))=p%y_face_fraction
      g%y_face_centroid_x(width+1:width+n(1),width:width+n(2),width+1:width+n(3))=p%y_face_centroid_x
      g%y_face_centroid_z(width+1:width+n(1),width:width+n(2),width+1:width+n(3))=p%y_face_centroid_z
      g%z_face_fraction(width+1:width+n(1),width+1:width+n(2),width:width+n(3))=p%z_face_fraction
      g%z_face_centroid_x(width+1:width+n(1),width+1:width+n(2),width:width+n(3))=p%z_face_centroid_x
      g%z_face_centroid_y(width+1:width+n(1),width+1:width+n(2),width:width+n(3))=p%z_face_centroid_y
      if(allocated(g%halo_volume_fraction)) g%halo_volume_fraction(1:g%nx,1:g%ny,1:g%nz)=g%volume_fraction
      end associate
      if(.not.candidate(node)%is_valid()) return
    end do
    call move_alloc(candidate,patches);ok=.true.
  end subroutine build_eb_tree_extended_geometry_3d

  subroutine initialize_reactive_eb_tree_3d(species,topology,finest,root_state,root_temperature,tree,ok,frame,internal_energy)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_patch_tree_3d), intent(in) :: topology
    type(eb_geometry_3d), intent(in) :: finest
    real(dp), intent(in) :: root_state(:,:,:,:),root_temperature(:,:,:)
    type(eb_reactive_tree_3d), intent(inout) :: tree
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy
    type(amr_tree_state_3d) :: held
    type(eb_reactive_tree_3d) :: candidate
    type(eb_geometry_3d), allocatable :: levels(:)
    real(dp), allocatable :: root_held(:,:,:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer :: count,level,node,parent,ratio,lo(3),hi(3),d(3),i,j,k,index(3),nv
    logical :: valid
    ok=.false.
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    call amr_tree_integer_layout_3d(topology,origin,domain,ratios,valid)
    if (.not.valid.or..not.finest%is_valid()) return
    count=size(ratios);nv=reactive_nvar(size(species))
    if (any(int([finest%nx,finest%ny,finest%nz],int64)/=domain(:,count))) return
    if (any(shape(root_state)/=[nv,topology%nodes(1)%shape])) return
    if (any(shape(root_temperature)/=topology%nodes(1)%shape)) return
    allocate(levels(count));levels(count)=finest
    do level=count,2,-1
      call coarsen_eb_geometry_3d(levels(level),ratios(level),levels(level-1),valid)
      if (.not.valid) return
    end do
    candidate%topology=topology
    allocate(candidate%nodes(size(topology%nodes)))
    if(present(internal_energy)) allocate(held%fields(size(topology%nodes)))
    do node=1,size(topology%nodes)
      level=topology%nodes(node)%level+1;d=topology%nodes(node)%shape
      lo=int(origin(:,node))+1;hi=lo+d-1
      call extract_eb_geometry_3d(levels(level),lo,hi,candidate%nodes(node)%geometry,valid)
      if (.not.valid) return
      allocate(candidate%nodes(node)%state(nv,d(1),d(2),d(3)),candidate%nodes(node)%temperature(d(1),d(2),d(3)))
      if(present(internal_energy)) allocate(held%fields(node)%values(1,d(1),d(2),d(3)),source=0._dp)
      if (node==1) then
        candidate%nodes(node)%state=root_state;candidate%nodes(node)%temperature=root_temperature
        if(present(internal_energy)) then
          ! Form UEINT from the supplied thermal seed before any total-minus-
          ! kinetic recovery. Its independent precision survives every child.
          allocate(root_held(d(1),d(2),d(3)),source=0._dp)
          call initialize_nasa7_internal_energy_patch(species,root_state,root_temperature, &
            candidate%nodes(1)%geometry%volume_fraction,root_held,valid)
          if(.not.valid) return
          call prolong_reactive_eb_state_pcm_3d(species,levels(1),1,root_state,root_temperature, &
            candidate%nodes(1)%state,candidate%nodes(1)%temperature,valid,frame,root_held, &
            held%fields(1)%values(1,:,:,:))
        else
          call prolong_reactive_eb_state_pcm_3d(species,levels(1),1,root_state,root_temperature, &
            candidate%nodes(1)%state,candidate%nodes(1)%temperature,valid,frame)
        endif
      else
        parent=topology%nodes(node)%parent;ratio=topology%nodes(node)%patch%refinement_ratio
        call bounds(topology,node,lo,hi)
        do k=1,d(3)
        do j=1,d(2)
        do i=1,d(1)
          index=lo+([i,j,k]-1)/ratio
          candidate%nodes(node)%state(:,i,j,k)=candidate%nodes(parent)%state(:,index(1),index(2),index(3))
          candidate%nodes(node)%temperature(i,j,k)=candidate%nodes(parent)%temperature(index(1),index(2),index(3))
        end do
        end do
        end do
        if(present(internal_energy)) then
          call prolong_reactive_eb_state_pcm_3d(species,candidate%nodes(node)%geometry,ratio, &
            candidate%nodes(parent)%state(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
            candidate%nodes(parent)%temperature(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
            candidate%nodes(node)%state,candidate%nodes(node)%temperature,valid,frame, &
            held%fields(parent)%values(1,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),held%fields(node)%values(1,:,:,:))
        else
          call prolong_reactive_eb_state_pcm_3d(species,candidate%nodes(node)%geometry,ratio, &
            candidate%nodes(parent)%state(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
            candidate%nodes(parent)%temperature(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
            candidate%nodes(node)%state,candidate%nodes(node)%temperature,valid,frame)
        endif
      end if
      if (.not.valid) return
    end do
    if (.not.valid_tree(candidate,nv)) return
    tree=candidate
    if(present(internal_energy)) call move_alloc(held%fields,internal_energy%fields)
    ok=.true.
  end subroutine initialize_reactive_eb_tree_3d

  ! Exposed fluid cells vote without averaging away fine temperature peaks.
  ! Reuse regular AMR's deterministic tile builder, then omit solid-only nodes.
  subroutine tag_reactive_eb_tree_temperature_3d(tree,finest,thresholds,ratios,tile_shape,buffer, &
      maximum_nodes,topology,ok)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: finest
    real(dp), intent(in) :: thresholds(:)
    integer, intent(in) :: ratios(:),tile_shape(3),buffer,maximum_nodes
    type(amr_patch_tree_3d), intent(inout) :: topology
    logical, intent(out) :: ok
    type(hot_cell), allocatable :: cells(:)
    type(amr_patch_tree_3d) :: candidate,pruned
    type(eb_geometry_3d), allocatable :: levels(:)
    type(eb_geometry_3d) :: section
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: old_ratios(:),shared(:),mapping(:)
    integer(int64) :: common(3),factor(3),total
    logical, allocatable :: active(:,:,:)
    integer :: nv,node,child,level,count,pass,i,j,k,n(3),lo(3),hi(3),used
    real(dp) :: temperature
    logical :: valid
    ok=.false.
    if (.not.allocated(tree%nodes)) return
    if (size(tree%nodes)<1) return
    if (.not.allocated(tree%nodes(1)%state)) return
    nv=size(tree%nodes(1)%state,1)
    if (.not.valid_tree(tree,nv).or..not.finest%is_valid()) return
    if (size(thresholds)/=size(ratios)) return
    if (.not.all(ieee_is_finite(thresholds)).or.any(thresholds<=0._dp)) return
    if (any(ratios<2).or.any(tile_shape<1).or.buffer<0.or.maximum_nodes<1) return
    do level=2,size(thresholds)
      if (thresholds(level)<thresholds(level-1)) return
    end do
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,old_ratios,valid)
    if (.not.valid) return
    level=min(size(old_ratios)-1,size(ratios))
    if (any(old_ratios(2:level+1)/=ratios(:level))) return
    count=max(size(old_ratios),size(ratios)+1)
    allocate(shared(count));shared(:size(old_ratios))=old_ratios
    shared(2:size(ratios)+1)=ratios;common=domain(:,1)
    do level=2,count
      if (any(common>(huge(1_int64)/4)/int(shared(level),int64))) return
      common=common*int(shared(level),int64)
    end do
    if (any(common/=int([finest%nx,finest%ny,finest%nz],int64))) return
    allocate(levels(count));levels(count)=finest
    do level=count,2,-1
      call coarsen_eb_geometry_3d(levels(level),shared(level),levels(level-1),valid)
      if (.not.valid) return
    end do
    do node=1,size(tree%nodes)
      level=tree%topology%nodes(node)%level+1
      lo=int(origin(:,node))+1;hi=lo+tree%topology%nodes(node)%shape-1
      call extract_eb_geometry_3d(levels(level),lo,hi,section,valid)
      if (.not.valid) return
      if (.not.matching_eb_metrics_3d(section,tree%nodes(node)%geometry)) return
    end do
    total=0
    do pass=1,2
      used=0
      do node=1,size(tree%nodes)
        n=tree%topology%nodes(node)%shape;level=tree%topology%nodes(node)%level+1
        factor=common/domain(:,level)
        allocate(active(n(1),n(2),n(3)));active=tree%nodes(node)%geometry%volume_fraction>0._dp
        do child=node+1,size(tree%nodes)
          if (tree%topology%nodes(child)%parent/=node) cycle
          call bounds(tree%topology,child,lo,hi)
          active(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
        end do
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if (.not.active(i,j,k)) cycle
          temperature=tree%nodes(node)%temperature(i,j,k)
          if (.not.ieee_is_finite(temperature).or.temperature<=0._dp) return
          if (size(thresholds)==0) cycle
          if (temperature<thresholds(1)) cycle
          if (pass==1) then
            total=total+1
            if (total>int(huge(1),int64)) return
          else
            used=used+1
            cells(used)%lower=(origin(:,node)+int([i,j,k]-1,int64))*factor
            cells(used)%upper=cells(used)%lower+factor;cells(used)%temperature=temperature
          end if
        end do
        end do
        end do
        deallocate(active)
      end do
      if (pass==1) allocate(cells(int(total)))
    end do
    call build_tagged_tree_3d(tree%topology,cells,common,thresholds,ratios,tile_shape,buffer, &
      maximum_nodes,[.false.,.false.,.false.],candidate,valid)
    if (.not.valid) return
    call amr_tree_integer_layout_3d(candidate,origin,domain,old_ratios,valid)
    if (.not.valid) return
    allocate(mapping(size(candidate%nodes)));mapping=0
    pruned=candidate;used=1;mapping(1)=1
    do node=2,size(candidate%nodes)
      child=candidate%nodes(node)%parent
      if (mapping(child)==0) cycle
      level=candidate%nodes(node)%level+1
      lo=int(origin(:,node))+1;hi=lo+candidate%nodes(node)%shape-1
      call extract_eb_geometry_3d(levels(level),lo,hi,section,valid)
      if (.not.valid) return
      if (.not.any(section%volume_fraction>0._dp)) cycle
      used=used+1;mapping(node)=used
      pruned%nodes(used)=candidate%nodes(node);pruned%nodes(used)%parent=mapping(child)
    end do
    pruned%nodes=pruned%nodes(:used)
    if (.not.pruned%is_valid()) return
    topology=pruned;ok=.true.
  end subroutine tag_reactive_eb_tree_temperature_3d

  ! Conservative remap on a shared refinement lattice. finest describes the
  ! deeper of old/new trees; ratios on shared levels and root cells must match.
  ! Averaged old cells supply coarsening, exact same-level overlap or ancestor
  ! PCM supplies new refinement. Geometry is immutable wherever old nodes exist.
  subroutine remap_reactive_eb_tree_3d(species,topology,finest,tree,ok,reaction_source,frame,heat_release,internal_energy)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_patch_tree_3d), intent(in) :: topology
    type(eb_geometry_3d), intent(in) :: finest
    type(eb_reactive_tree_3d), intent(inout) :: tree
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source,heat_release,internal_energy
    type(amr_tree_state_3d) :: old_source,new_source,old_heat,new_heat,old_held,new_held
    type(eb_reactive_tree_3d) :: snapshot,candidate,expected
    type(eb_geometry_3d), allocatable :: levels(:)
    type(eb_geometry_3d) :: section
    integer(int64), allocatable :: old_origin(:,:),old_domain(:,:),new_origin(:,:),new_domain(:,:)
    integer, allocatable :: old_ratios(:),new_ratios(:),ratios(:)
    integer(int64) :: global(3),mapped(3),local(3)
    integer :: nv,nlevels,old_levels,new_levels,level,node,donor,depth,i,j,k,d(3),index(3)
    real(dp) :: primitive(reactive_nprim(size(species))),temperature,speed
    logical :: valid,found
    ok=.false.;nv=reactive_nvar(size(species))
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    if (.not.valid_tree(tree,nv).or..not.finest%is_valid()) return
    if(present(reaction_source)) then
      if(.not.reaction_source%is_valid(tree%topology)) return
      if(size(reaction_source%fields(1)%values,1)/=nv-4) return
      old_source=reaction_source
      call average_down_reactive_eb_tree_field_3d(tree,old_source,valid)
      if(.not.valid) return
    end if
    if(present(heat_release)) then
      if(.not.heat_release%is_valid(tree%topology)) return
      if(size(heat_release%fields(1)%values,1)/=1) return
      old_heat=heat_release
      call average_down_reactive_eb_tree_field_3d(tree,old_heat,valid)
      if(.not.valid) return
    end if
    call amr_tree_integer_layout_3d(tree%topology,old_origin,old_domain,old_ratios,valid)
    if (.not.valid) return
    call amr_tree_integer_layout_3d(topology,new_origin,new_domain,new_ratios,valid)
    if (.not.valid) return
    old_levels=size(old_ratios);new_levels=size(new_ratios);nlevels=max(old_levels,new_levels)
    if (any(old_domain(:,1)/=new_domain(:,1))) return
    level=min(old_levels,new_levels)
    if (any(old_ratios(:level)/=new_ratios(:level))) return
    if (old_levels>=new_levels) then
      ratios=old_ratios
      if (any(int([finest%nx,finest%ny,finest%nz],int64)/=old_domain(:,nlevels))) return
    else
      ratios=new_ratios
      if (any(int([finest%nx,finest%ny,finest%nz],int64)/=new_domain(:,nlevels))) return
    end if
    allocate(levels(nlevels));levels(nlevels)=finest
    do level=nlevels,2,-1
      call coarsen_eb_geometry_3d(levels(level),ratios(level),levels(level-1),valid)
      if (.not.valid) return
    end do
    snapshot=tree
    if(present(internal_energy)) then
      old_held=internal_energy
      call average_down_reactive_eb_tree_3d(species,snapshot,valid,frame,old_held)
    else
      call average_down_reactive_eb_tree_3d(species,snapshot,valid,frame)
    endif
    if (.not.valid) return
    if(present(internal_energy)) then
      ! Build geometry without recovering root gas through unresolved E-K.
      ! Native gas/UEINT donors are copied together below, then restricted.
      do node=1,size(tree%nodes)
        level=tree%topology%nodes(node)%level+1;d=tree%topology%nodes(node)%shape
        call extract_eb_geometry_3d(levels(level),int(old_origin(:,node))+1, &
          int(old_origin(:,node))+d,section,valid)
        if(.not.valid) return
        if(.not.matching_eb_metrics_3d(tree%nodes(node)%geometry,section)) return
      enddo
      candidate%topology=topology
      allocate(candidate%nodes(size(topology%nodes)))
      do node=1,size(topology%nodes)
        level=topology%nodes(node)%level+1;d=topology%nodes(node)%shape
        call extract_eb_geometry_3d(levels(level),int(new_origin(:,node))+1, &
          int(new_origin(:,node))+d,candidate%nodes(node)%geometry,valid)
        if(.not.valid) return
        allocate(candidate%nodes(node)%state(nv,d(1),d(2),d(3)),candidate%nodes(node)%temperature(d(1),d(2),d(3)))
        candidate%nodes(node)%state=0;candidate%nodes(node)%temperature=0
      enddo
      allocate(new_held%fields(size(candidate%nodes)))
    else
    call initialize_reactive_eb_tree_3d(species,tree%topology,levels(old_levels),snapshot%nodes(1)%state, &
      snapshot%nodes(1)%temperature,expected,valid,frame)
    if (.not.valid) return
    do node=1,size(tree%nodes)
      if (.not.matching_eb_metrics_3d(tree%nodes(node)%geometry,expected%nodes(node)%geometry)) return
    end do
    call initialize_reactive_eb_tree_3d(species,topology,levels(new_levels),snapshot%nodes(1)%state, &
      snapshot%nodes(1)%temperature,candidate,valid,frame)
    if (.not.valid) return
    endif
    if(present(reaction_source)) allocate(new_source%fields(size(candidate%nodes)))
    if(present(heat_release)) allocate(new_heat%fields(size(candidate%nodes)))
    do node=1,size(candidate%nodes)
      if(present(reaction_source)) then
        d=topology%nodes(node)%shape
        allocate(new_source%fields(node)%values(nv-4,d(1),d(2),d(3)),source=0._dp)
      end if
      if(present(heat_release)) then
        d=topology%nodes(node)%shape
        allocate(new_heat%fields(node)%values(1,d(1),d(2),d(3)),source=0._dp)
      end if
      level=topology%nodes(node)%level+1;d=topology%nodes(node)%shape
      if(present(internal_energy)) allocate(new_held%fields(node)%values(1,d(1),d(2),d(3)),source=0._dp)
      do k=1,d(3)
      do j=1,d(2)
      do i=1,d(1)
        global=new_origin(:,node)+int([i,j,k]-1,int64);found=.false.
        do depth=min(level,old_levels),1,-1
          mapped=global/(new_domain(:,level)/old_domain(:,depth))
          do donor=1,size(snapshot%nodes)
            if (snapshot%topology%nodes(donor)%level+1/=depth) cycle
            local=mapped-old_origin(:,donor)+1_int64
            if (any(local<1_int64).or.any(local>int(snapshot%topology%nodes(donor)%shape,int64))) cycle
            index=int(local)
            if (candidate%nodes(node)%geometry%volume_fraction(i,j,k)>0._dp) then
              if (snapshot%nodes(donor)%geometry%volume_fraction(index(1),index(2),index(3))<=0._dp) return
              if(present(internal_energy)) then
              call recover_reactive_eb_cell_3d(species,candidate%nodes(node)%geometry,[i,j,k], &
                snapshot%nodes(donor)%state(:,index(1),index(2),index(3)), &
                snapshot%nodes(donor)%temperature(index(1),index(2),index(3)),primitive,temperature,speed,valid,frame, &
                old_held%fields(donor)%values(1,index(1),index(2),index(3)))
              else
              call recover_reactive_eb_cell_3d(species,candidate%nodes(node)%geometry,[i,j,k], &
                snapshot%nodes(donor)%state(:,index(1),index(2),index(3)), &
                snapshot%nodes(donor)%temperature(index(1),index(2),index(3)),primitive,temperature,speed,valid,frame)
              endif
              if (.not.valid) return
            end if
            candidate%nodes(node)%state(:,i,j,k)=snapshot%nodes(donor)%state(:,index(1),index(2),index(3))
            candidate%nodes(node)%temperature(i,j,k)=snapshot%nodes(donor)%temperature(index(1),index(2),index(3))
            if(present(internal_energy)) then
              new_held%fields(node)%values(1,i,j,k)=old_held%fields(donor)%values(1,index(1),index(2),index(3))
              if(candidate%nodes(node)%geometry%volume_fraction(i,j,k)>0) &
                candidate%nodes(node)%temperature(i,j,k)=temperature
            endif
            ! A copied rotating energy has a different thermal part at a new
            ! cell center. Same-level overlap retains the exact stored value.
            if(present(frame).and.depth/=level) then
              if(candidate%nodes(node)%geometry%volume_fraction(i,j,k)>0) &
                candidate%nodes(node)%temperature(i,j,k)=temperature
            end if
            if(present(reaction_source)) then
              if(candidate%nodes(node)%geometry%volume_fraction(i,j,k)>0) &
                new_source%fields(node)%values(:,i,j,k)=old_source%fields(donor)%values(:,index(1),index(2),index(3))
            end if
            if(present(heat_release)) then
              if(candidate%nodes(node)%geometry%volume_fraction(i,j,k)>0) &
                new_heat%fields(node)%values(:,i,j,k)=old_heat%fields(donor)%values(:,index(1),index(2),index(3))
            end if
            found=.true.;exit
          end do
          if (found) exit
        end do
        if (.not.found) return
      end do
      end do
      end do
    end do
    if(present(internal_energy)) then
      call average_down_reactive_eb_tree_3d(species,candidate,valid,frame,new_held)
    else
      call average_down_reactive_eb_tree_3d(species,candidate,valid,frame)
    endif
    if (.not.valid) return
    if(present(reaction_source)) then
      call average_down_reactive_eb_tree_field_3d(candidate,new_source,valid)
      if(.not.valid) return
    end if
    if(present(heat_release)) then
      call average_down_reactive_eb_tree_field_3d(candidate,new_heat,valid)
      if(.not.valid) return
    end if
    if(present(reaction_source)) reaction_source=new_source
    if(present(heat_release)) heat_release=new_heat
    tree=candidate
    if(present(internal_energy)) call move_alloc(new_held%fields,internal_energy%fields)
    ok=.true.
  end subroutine remap_reactive_eb_tree_3d

  ! Signed cell sources follow the same EB-volume average as the pinned
  ! PeleC Reactions_Type. Publish only after every descendant succeeds.
  subroutine average_down_reactive_eb_tree_field_3d(layout,field,good,minimum_depth)
    type(eb_reactive_tree_3d), intent(in) :: layout
    type(amr_tree_state_3d), intent(inout) :: field
    logical, intent(out) :: good
    integer, optional, intent(in) :: minimum_depth
    type(amr_tree_state_3d) :: candidate
    integer :: child,parent,r,lo(3),hi(3),x,y,z,first(3),last(3),component,nv,p,depth
    real(dp) :: volume
    real(dp), allocatable :: value(:)
    good=.false.
    if(.not.field%is_valid(layout%topology)) return
    if(.not.allocated(layout%nodes)) return
    if(size(layout%nodes)/=size(layout%topology%nodes)) return
    do p=1,size(layout%nodes)
      if(.not.layout%nodes(p)%geometry%is_valid()) return
      if(any([layout%nodes(p)%geometry%nx,layout%nodes(p)%geometry%ny, &
        layout%nodes(p)%geometry%nz]/=layout%topology%nodes(p)%shape)) return
    enddo
    depth=0
    if(present(minimum_depth)) depth=minimum_depth
    if(depth<0.or.depth>=layout%topology%level_count()) return
    nv=size(field%fields(1)%values,1)
    allocate(value(nv));candidate=field
    do child=size(layout%nodes),2,-1
      parent=layout%topology%nodes(child)%parent
      if(layout%topology%nodes(parent)%level<depth) cycle
      r=layout%topology%nodes(child)%patch%refinement_ratio
      call bounds(layout%topology,child,lo,hi)
      do z=lo(3),hi(3)
      do y=lo(2),hi(2)
      do x=lo(1),hi(1)
        first=([x,y,z]-lo)*r+1;last=first+r-1
        associate(v=>layout%nodes(child)%geometry%volume_fraction(first(1):last(1),first(2):last(2),first(3):last(3)))
        volume=sum(v);value=0
        if(volume>0) then
          do component=1,nv
            value(component)=sum(v*candidate%fields(child)%values(component, &
              first(1):last(1),first(2):last(2),first(3):last(3)))/volume
          end do
        end if
        end associate
        if(.not.all(ieee_is_finite(value))) return
        candidate%fields(parent)%values(:,x,y,z)=value
      end do
      end do
      end do
    end do
    field=candidate;good=.true.
  end subroutine average_down_reactive_eb_tree_field_3d


  logical function valid_tree(tree,nv) result(valid)
    type(eb_reactive_tree_3d), intent(in) :: tree
    integer, intent(in) :: nv
    type(eb_geometry_3d) :: section,restricted
    integer :: node,parent,lo(3),hi(3),d(3)
    logical :: good
    valid=.false.
    if (.not.tree%topology%is_valid().or..not.allocated(tree%nodes)) return
    if (size(tree%nodes)/=size(tree%topology%nodes)) return
    do node=1,size(tree%nodes)
      if (.not.tree%nodes(node)%geometry%is_valid()) return
      if (.not.allocated(tree%nodes(node)%state).or..not.allocated(tree%nodes(node)%temperature)) return
      d=tree%topology%nodes(node)%shape
      if (any([tree%nodes(node)%geometry%nx,tree%nodes(node)%geometry%ny,tree%nodes(node)%geometry%nz]/=d)) return
      if (any(shape(tree%nodes(node)%state)/=[nv,d]).or.any(shape(tree%nodes(node)%temperature)/=d)) return
      if (node==1) cycle
      parent=tree%topology%nodes(node)%parent
      call bounds(tree%topology,node,lo,hi)
      call coarsen_eb_geometry_3d(tree%nodes(node)%geometry,tree%topology%nodes(node)%patch%refinement_ratio,restricted,good)
      if (.not.good) return
      call extract_eb_geometry_3d(tree%nodes(parent)%geometry,lo,hi,section,good)
      if (.not.good) return
      if (.not.matching_eb_metrics_3d(section,restricted)) return
    end do
    valid=.true.
  end function valid_tree

  subroutine average_down_reactive_eb_tree_3d(species,tree,ok,frame,internal_energy,direct_children)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(eb_reactive_tree_3d), intent(inout) :: tree
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy
    logical, optional, intent(in) :: direct_children
    type(amr_tree_state_3d) :: held
    type(eb_reactive_tree_3d) :: candidate
    integer :: node,parent,lo(3),hi(3)
    logical :: valid
    ok=.false.
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    if (.not.valid_tree(tree,reactive_nvar(size(species)))) return
    if(present(internal_energy)) then
      if(.not.internal_energy%is_valid(tree%topology)) return
      do node=1,size(tree%nodes)
        if(size(internal_energy%fields(node)%values,1)/=1) return
        if(any(.not.ieee_is_finite(internal_energy%fields(node)%values(1,:,:,:)).and. &
          tree%nodes(node)%geometry%volume_fraction>0)) return
      enddo
      held=internal_energy
    endif
    candidate=tree
    do node=size(tree%nodes),2,-1
      parent=tree%topology%nodes(node)%parent
      call bounds(tree%topology,node,lo,hi)
      if(present(direct_children)) then
        if(direct_children) then
      if(present(internal_energy)) then
      call restrict_reactive_eb_state_3d(species,tree%nodes(node)%geometry, &
        tree%topology%nodes(node)%patch%refinement_ratio,tree%nodes(node)%state,tree%nodes(node)%temperature, &
        candidate%nodes(parent)%state(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
        candidate%nodes(parent)%temperature(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),valid,frame, &
        internal_energy%fields(node)%values(1,:,:,:),held%fields(parent)%values(1,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      else
      call restrict_reactive_eb_state_3d(species,tree%nodes(node)%geometry, &
        tree%topology%nodes(node)%patch%refinement_ratio,tree%nodes(node)%state,tree%nodes(node)%temperature, &
        candidate%nodes(parent)%state(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
        candidate%nodes(parent)%temperature(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),valid,frame)
      endif
          if(.not.valid) return
          cycle
        endif
      endif
      if(present(internal_energy)) then
      call restrict_reactive_eb_state_3d(species,candidate%nodes(node)%geometry, &
        tree%topology%nodes(node)%patch%refinement_ratio,candidate%nodes(node)%state,candidate%nodes(node)%temperature, &
        candidate%nodes(parent)%state(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
        candidate%nodes(parent)%temperature(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),valid,frame, &
        held%fields(node)%values(1,:,:,:),held%fields(parent)%values(1,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)))
      else
      call restrict_reactive_eb_state_3d(species,candidate%nodes(node)%geometry, &
        tree%topology%nodes(node)%patch%refinement_ratio,candidate%nodes(node)%state,candidate%nodes(node)%temperature, &
        candidate%nodes(parent)%state(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)), &
        candidate%nodes(parent)%temperature(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),valid,frame)
      endif
      if (.not.valid) return
    end do
    tree=candidate
    if(present(internal_energy)) call move_alloc(held%fields,internal_energy%fields)
    ok=.true.
  end subroutine average_down_reactive_eb_tree_3d

  ! Bound every patch update, including hidden parent cells which are advanced
  ! before restriction. Global donor maxima cover peer and ancestor halos.
  subroutine reactive_eb_tree_timestep_3d(species,transport,tree,hydro_cfl,transport_cfl, &
      transport_enabled,viscosity,conduction,diffusion,dt,ok,subcycling, &
      inflow_mach,inflow_temperature,inflow_pressure,boundary_faces,inflow_state,wall_temperature,coupled_subcycling,frame, &
      constant_coefficients,internal_energy)
    real(dp), optional, intent(in) :: constant_coefficients(:)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy
    type(amr_tree_state_3d) :: held
    real(dp), intent(in) :: hydro_cfl,transport_cfl
    logical, intent(in) :: transport_enabled,viscosity,conduction,diffusion
    real(dp), intent(inout) :: dt
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: wall_temperature(6)
    character(len=16) :: thermal_faces(6)
    integer :: direction
    logical :: thermal_wall
    logical, intent(out) :: ok
    logical, optional, intent(in) :: subcycling,coupled_subcycling
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    type(eb_reactive_tree_3d) :: snapshot
    real(dp) :: speeds(4),donors(4),diffusivity,maximum_diffusivity,interval,candidate,factor
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:)
    real(dp), allocatable :: ordinary(:,:,:,:)
    real(dp), allocatable :: patch_held(:,:,:),ghost_held(:,:,:)
    integer :: node,level,n(3)
    logical :: valid,active,recursive,coupled
    ok=.false.
    if(present(constant_coefficients)) then
      if(.not.valid_constant_transport_coefficients(constant_coefficients)) return
    endif
    if (.not.all(ieee_is_finite([hydro_cfl,transport_cfl]))) return
    if (hydro_cfl<=0._dp.or.hydro_cfl>1._dp.or.transport_cfl<=0._dp.or.transport_cfl>.5_dp) return
    if (.not.valid_tree(tree,reactive_nvar(size(species)))) return
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    thermal_wall=.false.
    if(present(wall_temperature)) then
      if(.not.present(boundary_faces)) return
      if(.not.all(ieee_is_finite(wall_temperature)).or.any(wall_temperature<=0)) return
    end if
    if(present(boundary_faces)) thermal_wall=any(boundary_faces=='isothermal_wall')
    if(thermal_wall) then
      if(.not.present(wall_temperature).or..not.transport_enabled.or..not.conduction) return
    end if
    coupled=.false.
    if(present(coupled_subcycling)) coupled=coupled_subcycling
    recursive=.false.
    if (present(subcycling)) recursive=subcycling
    if (recursive.or.thermal_wall) then
      call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
      if (.not.valid) return
    end if
    snapshot=tree
    if(present(internal_energy)) then
      held=internal_energy
      call average_down_reactive_eb_tree_3d(species,snapshot,valid,frame,held)
    else
      call average_down_reactive_eb_tree_3d(species,snapshot,valid,frame)
    endif
    if (.not.valid) return
    donors=0._dp;maximum_diffusivity=0._dp;candidate=huge(1._dp);active=.false.
    do node=1,size(tree%nodes)
      associate(g=>snapshot%nodes(node)%geometry,u=>snapshot%nodes(node)%state,t=>snapshot%nodes(node)%temperature)
      if (.not.any(g%volume_fraction>0._dp)) cycle
      active=.true.
      ordinary=u
      if(present(internal_energy)) then
        patch_held=held%fields(node)%values(1,:,:,:)
      else
        call ordinary_transport_energy(frame,g,[1,1,1],g%volume_fraction>0,ordinary,valid)
      endif
      if(.not.valid) return
      call reactive_eb_signal_speeds_3d(species,ordinary,t,g,speeds,valid,patch_held)
      if (.not.valid) return
      donors=max(donors,speeds)
      if (present(inflow_mach).or.present(inflow_temperature).or.present(inflow_pressure).or. &
        present(boundary_faces).or.present(inflow_state)) then
        n=tree%topology%nodes(node)%shape
        allocate(ghost(size(u,1),0:n(1)+1,0:n(2)+1,0:n(3)+1),gt(0:n(1)+1,0:n(2)+1,0:n(3)+1))
        if(present(internal_energy)) then
          allocate(ghost_held(0:n(1)+1,0:n(2)+1,0:n(3)+1))
          call fill_reactive_eb_tree_halo_3d(species,snapshot,node,ghost,gt,valid, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
            boundary_faces=boundary_faces,inflow_state=inflow_state,frame=frame,restrict_snapshot=.false., &
            internal_energy=held,ghost_internal_energy=ghost_held)
        else
          call fill_reactive_eb_tree_halo_3d(species,snapshot,node,ghost,gt,valid, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
            boundary_faces=boundary_faces,inflow_state=inflow_state,frame=frame)
        endif
        if (.not.valid) return
        if(.not.present(internal_energy)) call ordinary_transport_energy(frame,g,[0,0,0],ghost(1,:,:,:)>0,ghost,valid)
        if(.not.valid) return
        call reactive_eb_halo_bounds_3d(species,transport,g,ghost,gt,transport_enabled, &
          viscosity,conduction,diffusion,speeds,diffusivity,valid, &
          constant_coefficients=constant_coefficients,internal_energy_density=ghost_held)
        if (.not.valid) return
        donors=max(donors,speeds);maximum_diffusivity=max(maximum_diffusivity,diffusivity)
        deallocate(ghost,gt)
        if(allocated(ghost_held)) deallocate(ghost_held)
      end if
      if (transport_enabled) then
        call reactive_eb_transport_patch_timestep_3d(species,transport,ordinary,t,g,transport_cfl, &
          viscosity,conduction,diffusion,interval,diffusivity,valid, &
          constant_coefficients=constant_coefficients,internal_energy_density=patch_held)
        if (.not.valid) return
        maximum_diffusivity=max(maximum_diffusivity,diffusivity)
      end if
      end associate
    end do
    if (.not.active) return
    do node=1,size(tree%nodes)
      associate(g=>snapshot%nodes(node)%geometry,u=>snapshot%nodes(node)%state,t=>snapshot%nodes(node)%temperature)
      if (.not.any(g%volume_fraction>0._dp)) cycle
      ordinary=u
      if(present(internal_energy)) then
        patch_held=held%fields(node)%values(1,:,:,:)
      else
        call ordinary_transport_energy(frame,g,[1,1,1],g%volume_fraction>0,ordinary,valid)
      endif
      if(.not.valid) return
      call compute_reactive_eb_euler_cfl_timestep_3d(species,ordinary,t,g,hydro_cfl,interval,valid,donors,patch_held)
      if (.not.valid) return
      factor=1._dp
      if (recursive) then
        level=tree%topology%nodes(node)%level+1
        factor=real(domain(1,level),dp)/real(domain(1,1),dp)
      end if
      ! Compare before multiplying so an inactive huge bound cannot overflow.
      if (interval<candidate/factor) candidate=interval*factor
      if (transport_enabled) then
        call reactive_eb_transport_patch_timestep_3d(species,transport,ordinary,t,g,transport_cfl, &
          viscosity,conduction,diffusion,interval,diffusivity,valid,maximum_diffusivity, &
          constant_coefficients=constant_coefficients,internal_energy_density=patch_held)
        if (.not.valid) return
        ! Conservatively use a full-interval diffusion bound for split halves.
        if(.not.coupled) factor=factor*factor
        if (interval<candidate/factor) candidate=interval*factor
        if(thermal_wall) then
          level=tree%topology%nodes(node)%level+1;n=tree%topology%nodes(node)%shape
          thermal_faces='outflow'
          do direction=1,3
            if(origin(direction,node)==0) thermal_faces(2*direction-1)=boundary_faces(2*direction-1)
            if(origin(direction,node)+n(direction)==domain(direction,level)) &
              thermal_faces(2*direction)=boundary_faces(2*direction)
          end do
          call reactive_eb_isothermal_timestep_3d(species,transport,ordinary,t,g,thermal_faces, &
            wall_temperature,transport_cfl,interval,valid, &
            constant_coefficients=constant_coefficients,internal_energy_density=patch_held)
          if(.not.valid) return
          if(interval<candidate/factor) candidate=interval*factor
        end if
      end if
      end associate
    end do
    dt=candidate;ok=.true.
  end subroutine reactive_eb_tree_timestep_3d

  subroutine reset_held_eb_tree_family_3d(species,gas,stored,delta,policy,depth,change,good,frame, &
    leaf_only,internal_energy_change,reset_level,covered_only)
    ! -1 (default): reset normally; -2: validate on entry, repairing only
    ! failed held EOS states; >=0: reset that level and validate the others.
    ! Ordinary reacting MOL retains the fine-level total/held difference
    ! left by parent reflux until that level advances. A valid held EOS does
    ! not authorize overwriting conserved energy at another clock's boundary.
    integer, optional, intent(in) :: reset_level
    logical :: validate_first
    type(nasa7_species), intent(in) :: species(:)
    type(dual_energy_policy), intent(in) :: policy
    integer, intent(in) :: depth
    type(rotating_frame_3d), optional, intent(in) :: frame
    logical, optional, intent(in) :: leaf_only,covered_only
    type(eb_reactive_tree_3d), intent(inout) :: gas
    type(amr_tree_state_3d), intent(inout) :: stored
    type(amr_tree_state_3d), intent(out) :: delta
    real(dp), intent(out) :: change
    real(dp), optional, intent(inout) :: internal_energy_change
    logical, intent(out) :: good
    real(dp) :: cell_change,phi,position(3),volume,pv(reactive_nprim(size(species))),temp,sound
    real(dp) :: held_before,held_cell_change,local_held_change
    integer :: p,d(3),x,y,z
    good=.false.;change=0;local_held_change=0
    allocate(delta%fields(size(gas%nodes)))
    do p=1,size(gas%nodes)
      allocate(delta%fields(p)%values,source=gas%nodes(p)%state)
      delta%fields(p)%values=0
      if(depth>=0.and.gas%topology%nodes(p)%level/=depth) cycle
      d=gas%topology%nodes(p)%shape
      associate(g=>gas%nodes(p)%geometry)
      volume=g%dx*g%dy*g%dz
      do z=1,d(3)
      do y=1,d(2)
      do x=1,d(1)
        if(g%volume_fraction(x,y,z)<=0) cycle
        if(present(leaf_only)) then
          if(leaf_only.and.covered(gas%topology,p,[x,y,z])) cycle
        endif
        if(present(covered_only)) then
          if(covered_only.and..not.covered(gas%topology,p,[x,y,z])) cycle
        endif
        phi=0
        if(present(frame)) then
          position=[g%x_lower,g%y_lower,g%z_lower]+(real([x,y,z],dp)-.5_dp)*[g%dx,g%dy,g%dz]
          call rotating_potential_3d(frame,position,phi,good)
          if(.not.good) return
        endif
        validate_first=.false.
        if(present(reset_level)) validate_first=reset_level==-2.or. &
          (reset_level>=0.and.gas%topology%nodes(p)%level/=reset_level)
        if(validate_first) then
          call recover_reactive_eb_cell_3d(species,g,[x,y,z],gas%nodes(p)%state(:,x,y,z), &
            gas%nodes(p)%temperature(x,y,z),pv,temp,sound,good,frame,stored%fields(p)%values(1,x,y,z))
          if(good) then
            ! Admission preserves the native clock's inverse-EOS seed.
            ! With the pinned stopping rule, recomputing a valid seed changes
            ! later pressure even when conserved and held energy are unchanged.
            if(reset_level/=-2) gas%nodes(p)%temperature(x,y,z)=temp
            cycle
          endif
          ! Fall through to the existing policy repair and ledger accounting.
        endif
        cell_change=0
        held_before=0
        if(present(internal_energy_change)) held_before=stored%fields(p)%values(1,x,y,z)
        call reset_nasa7_dual_energy_cell(species,gas%nodes(p)%state(:,x,y,z), &
          stored%fields(p)%values(1,x,y,z),gas%nodes(p)%temperature(x,y,z),policy,cell_change,good,phi)
        if(.not.good) return
        if(present(internal_energy_change)) then
          call compute_reactive_eb_held_projection(held_before,held_before,held_before, &
            stored%fields(p)%values(1,x,y,z),held_cell_change,good)
          if(.not.good) return
          if(.not.covered(gas%topology,p,[x,y,z])) &
            local_held_change=local_held_change+volume*g%volume_fraction(x,y,z)*held_cell_change
        endif
        delta%fields(p)%values(5,x,y,z)=cell_change
        if(.not.covered(gas%topology,p,[x,y,z])) change=change+volume*g%volume_fraction(x,y,z)*cell_change
        call recover_reactive_eb_cell_3d(species,g,[x,y,z],gas%nodes(p)%state(:,x,y,z), &
          gas%nodes(p)%temperature(x,y,z),pv,temp,sound,good,frame,stored%fields(p)%values(1,x,y,z))
        if(.not.good) return
        gas%nodes(p)%temperature(x,y,z)=temp
      enddo
      enddo
      enddo
      end associate
    enddo
    good=ieee_is_finite(change).and.ieee_is_finite(local_held_change)
    if(good.and.present(internal_energy_change)) internal_energy_change=local_held_change
  end subroutine reset_held_eb_tree_family_3d

  ! Close the matching accepted gas and UEINT families privately, before
  ! temperature admission. Only exposed leaves charge the numerical reset;
  ! their covered ancestors are rebuilt by average-down without a second charge.
  subroutine close_held_prepared_eb_tree_3d(species,tree,stored,flux,scalar_flux,policy,prefix,depth,ok, &
      boundary_faces,frame,extended_geometry,registers,scalar_registers,internal_energy_change,reset_current_level, &
      diagnostic_raw_increment,diagnostic_applied_increment)
    type(amr_tree_state_3d), optional, intent(inout) :: diagnostic_raw_increment,diagnostic_applied_increment
    type(amr_tree_state_3d) :: captured_raw
    logical, optional, intent(in) :: reset_current_level
    integer :: selected_reset_level
    type(nasa7_species), intent(in) :: species(:)
    type(eb_reactive_tree_3d), intent(inout) :: tree
    type(amr_tree_state_3d), intent(inout) :: stored
    type(eb_integrated_flux_3d), intent(inout) :: flux(:),scalar_flux(:)
    type(dual_energy_policy), intent(in) :: policy
    real(dp), intent(inout) :: prefix
    real(dp), optional, intent(inout) :: internal_energy_change
    integer, intent(in) :: depth
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry(:)
    type(state_redist_register_3d), optional, intent(inout) :: registers(:),scalar_registers(:)
    type(eb_reactive_tree_3d) :: candidate
    type(amr_tree_state_3d) :: held,correction,spill,increment,reset_delta,scalar_correction,scalar_spill
    type(eb_integrated_flux_3d), allocatable :: working(:),scalar_work(:),repair_faces(:),scalar_input(:)
    type(state_redist_register_3d), allocatable :: held_register_work(:)
    type(amr_tree_state_3d) :: held_before,joint_correction,joint_spill
    type(state_redist_register_3d), allocatable :: flux_registers(:)
    real(dp), allocatable :: joint_source(:,:,:,:)
    integer :: h
    logical :: flux_rereflux
    real(dp) :: change,new_prefix
    real(dp) :: held_reset_change,held_prefix
    integer :: nv,p,n(3),x,y,z,first_depth,last_depth,closure_attempt,scalar_first_depth
    logical :: good,has_spill,has_scalar_spill,closed,cell_ok
    ok=.false.;nv=reactive_nvar(size(species))
    selected_reset_level=-1
    if(present(reset_current_level)) then
      if(reset_current_level) selected_reset_level=depth
    endif
    if(.not.ieee_is_finite(prefix)) return
    if(present(registers).neqv.present(scalar_registers)) return
    if(present(registers).and..not.present(extended_geometry)) return
    flux_rereflux=present(extended_geometry).and..not.present(registers).and.depth+1<tree%topology%level_count()
    if(flux_rereflux) allocate(flux_registers(size(tree%nodes)))
    if(.not.stored%is_valid(tree%topology)) return
    if(size(stored%fields(1)%values,1)/=1) return
    candidate=tree;held=stored;repair_faces=flux;scalar_input=scalar_flux
    if(present(internal_energy_change)) then
      if(.not.ieee_is_finite(internal_energy_change)) return
    endif
    held_reset_change=0
    do closure_attempt=1,2
    scalar_work=scalar_input
    held_before=held
    first_depth=depth;last_depth=min(depth+1,tree%topology%level_count()-1)
    if(present(registers)) then
      call assemble_state_redist_interface_rereflux_3d(tree%topology,extended_geometry,registers,nv,depth, &
        correction,good,boundary_faces,spill,has_spill)
      if(.not.good) return
      if(has_spill) then
        first_depth=depth-1
        do p=1,size(tree%nodes)
          if(tree%topology%nodes(p)%level==first_depth) correction%fields(p)=spill%fields(p)
        enddo
      endif
      call assemble_state_redist_interface_rereflux_3d(tree%topology,extended_geometry,scalar_registers,1,depth, &
        scalar_correction,good,boundary_faces,scalar_spill,has_scalar_spill)
      if(.not.good) return
      scalar_first_depth=depth
      if(has_scalar_spill) then
        scalar_first_depth=depth-1
        do p=1,size(tree%nodes)
          if(tree%topology%nodes(p)%level==scalar_first_depth) scalar_correction%fields(p)=scalar_spill%fields(p)
        enddo
      endif
      first_depth=min(first_depth,scalar_first_depth)
      call assemble_reactive_eb_tree_reflux_3d(species,candidate,repair_faces,increment,working,good, &
        depth,last_depth,boundary_faces,frame,correction,first_depth)
      if(.not.good) return
      held_register_work=scalar_registers
      call synchronize_reactive_eb_tree_scalar_register_flux_3d(species,candidate,held,scalar_work, &
        extended_geometry,held_register_work,depth,good,boundary_faces)
    else
      call assemble_reactive_eb_tree_reflux_3d(species,candidate,repair_faces,increment,working,good, &
        depth,last_depth,boundary_faces,frame)
      if(.not.good) return
      call synchronize_reactive_eb_tree_scalar_flux_3d(species,candidate,held,scalar_work,good, &
        depth,last_depth,boundary_faces)
    endif
    if(.not.good) return
    if(present(diagnostic_raw_increment)) captured_raw=increment
    if(flux_rereflux) then
      ! Gas and actual UEINT share one spatial redistribution. The raw
      ! synchronization supplies the conservative mismatch, not a new source.
      do p=1,size(tree%nodes)
        n=tree%topology%nodes(p)%shape
        allocate(joint_source(nv+1,n(1),n(2),n(3)))
        joint_source(:nv,:,:,:)=increment%fields(p)%values
        joint_source(nv+1:nv+1,:,:,:)=held%fields(p)%values-held_before%fields(p)%values
        call initialize_flux_rereflux_register_3d(tree%topology,p,extended_geometry(p),joint_source, &
          flux_registers(p),good,boundary_faces)
        if(.not.good) return
        deallocate(joint_source)
      enddo
      call assemble_state_redist_interface_rereflux_3d(tree%topology,extended_geometry,flux_registers,nv+1,depth, &
        joint_correction,good,boundary_faces,joint_spill,has_spill)
      if(.not.good) return
      first_depth=depth
      if(has_spill) first_depth=depth-1
      correction=increment;scalar_correction=held_before
      do p=1,size(tree%nodes)
        n=tree%topology%nodes(p)%shape
        if(tree%topology%nodes(p)%level==depth) then
          h=(extended_geometry(p)%nx-n(1))/2
          joint_correction%fields(p)%values=joint_correction%fields(p)%values- &
            flux_registers(p)%coarse_delta(:,h+1:h+n(1),h+1:h+n(2),h+1:h+n(3))
        else if(has_spill.and.tree%topology%nodes(p)%level==first_depth) then
          joint_correction%fields(p)=joint_spill%fields(p)
        endif
        correction%fields(p)%values=joint_correction%fields(p)%values(:nv,:,:,:)
        scalar_correction%fields(p)%values=joint_correction%fields(p)%values(nv+1:nv+1,:,:,:)
      enddo
      ! Replace the raw mismatch by its re-refluxed form through the same
      ! routing used by spatial registers, including fine-covered recipients.
      call assemble_reactive_eb_tree_reflux_3d(species,candidate,repair_faces,increment,working,good, &
        depth,last_depth,boundary_faces,frame,correction,first_depth)
      if(.not.good) return
      held=held_before;scalar_work=scalar_input
      call synchronize_reactive_eb_tree_scalar_flux_3d(species,candidate,held,scalar_work,good, &
        depth,last_depth,boundary_faces,scalar_correction,first_depth)
      if(.not.good) return
    endif
    closed=.true.
    do p=1,size(tree%nodes)
      n=tree%topology%nodes(p)%shape
      do z=1,n(3)
      do y=1,n(2)
      do x=1,n(1)
        if(tree%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
        call add_eb_tree_delta_cell_3d(candidate%nodes(p)%state(:,x,y,z), &
          increment%fields(p)%values(:,x,y,z),cell_ok)
        if(.not.cell_ok) return
        closed=closed.and.all(candidate%nodes(p)%state(6:,x,y,z)>=0)
      enddo
      enddo
      enddo
    enddo
    if(closed) then
      call reset_held_eb_tree_family_3d(species,candidate,held,reset_delta,policy,-1,change,good,frame,leaf_only=.true., &
        internal_energy_change=held_reset_change,reset_level=selected_reset_level)
      closed=good
    endif
    if(closed) exit
    if(closure_attempt/=1) return
    candidate=tree;held=stored;repair_faces=flux;scalar_input=scalar_flux
    if(present(registers).or.flux_rereflux) then
      call repair_reactive_eb_tree_held_interface_budget_3d(species,candidate,held,repair_faces,scalar_input, &
        good,depth,boundary_faces,frame,correction,scalar_correction,rereflux_first_level=first_depth)
    else
    call repair_reactive_eb_tree_held_interface_budget_3d(species,candidate,held,repair_faces,scalar_input, &
      good,depth,boundary_faces,frame)
    endif
    if(.not.good) return
    enddo
    if(.not.closed) return
    held_prefix=held_reset_change
    if(present(internal_energy_change)) held_prefix=internal_energy_change+held_reset_change
    if(.not.ieee_is_finite(held_prefix)) return
    new_prefix=prefix+change
    if(.not.ieee_is_finite(new_prefix)) return
    call average_down_reactive_eb_tree_3d(species,candidate,good,frame,held)
    if(.not.good) return
    if(selected_reset_level>=0) then
      ! Covered parent cells supply future interpolation but carry no
      ! composite inventory. Repair them only after the closing average-down.
      call reset_held_eb_tree_family_3d(species,candidate,held,reset_delta,policy,selected_reset_level, &
        change,good,frame,covered_only=.true.)
      if(.not.good.or.change/=0) return
    endif
    tree=candidate;call move_alloc(held%fields,stored%fields)
    flux=working;scalar_flux=scalar_work;prefix=new_prefix
    if(present(registers)) then
      scalar_registers=held_register_work
      do p=1,size(tree%nodes)
        if(tree%topology%nodes(p)%level==depth) then
          if(allocated(registers(p)%coarse_delta)) registers(p)%coarse_delta=0
        else if(tree%topology%nodes(p)%level==depth+1) then
          if(allocated(registers(p)%fine_delta)) registers(p)%fine_delta=0
        endif
      enddo
    endif
    if(present(internal_energy_change)) internal_energy_change=held_prefix
    if(present(diagnostic_raw_increment)) diagnostic_raw_increment=captured_raw
    if(present(diagnostic_applied_increment)) diagnostic_applied_increment=increment
    ok=.true.
  end subroutine close_held_prepared_eb_tree_3d

  ! Reconstruct gas reaction forcing and the independent held spatial family.
  ! Density/species volume sources are available directly; avoid cancellation
  ! of resolved trace rates by subtracting two almost equal mass fluxes.
  subroutine reconstruct_reactive_eb_mol_held_interval_3d(species,tree,history,flux,scalar_flux,rates,scalar_rates,ok,frame)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(eb_mol_held_interval_3d), intent(in) :: history
    type(eb_integrated_flux_3d), intent(in) :: flux(:),scalar_flux(:)
    type(amr_tree_state_3d), intent(inout) :: rates,scalar_rates
    logical, intent(out) :: ok
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d) :: work,held_work
    real(dp) :: change(reactive_nvar(size(species))),new_div(reactive_nvar(size(species)))
    real(dp) :: held_change(1),volume
    integer :: nv,p,x,y,z
    ok=.false.;nv=reactive_nvar(size(species))
    if(history%prepared_redistribution) return
    if(.not.valid_held_mol_record_3d(tree,history,nv,frame)) return
    if(size(flux)/=size(tree%nodes).or.size(scalar_flux)/=size(tree%nodes)) return
    work=history%forcing;held_work=history%internal_energy_spatial_rate
    do p=1,size(tree%nodes)
      if(tree%topology%nodes(p)%level/=history%depth) then
        if(history%has_source_funding) then
          if(.not.valid_eb_integrated_flux_3d(flux(p),tree%nodes(p)%geometry,nv).or. &
            .not.valid_eb_integrated_flux_3d(scalar_flux(p),tree%nodes(p)%geometry,1)) return
        endif
        work%fields(p)%values=0;held_work%fields(p)%values=0;cycle
      endif
      associate(g=>tree%nodes(p)%geometry)
      if(.not.valid_eb_integrated_flux_3d(flux(p),g,nv).or. &
        .not.valid_eb_integrated_flux_3d(scalar_flux(p),g,1)) return
      do z=1,g%nz
      do y=1,g%ny
      do x=1,g%nx
        if(g%volume_fraction(x,y,z)<=0) then
          work%fields(p)%values(:,x,y,z)=0;held_work%fields(p)%values(:,x,y,z)=0;cycle
        endif
        volume=g%dx*g%dy*g%dz*g%volume_fraction(x,y,z)
        new_div=held_mol_face_impulse_3d(g,flux(p),x,y,z)
        change=held_mol_face_change_3d(g,flux(p),history%spatial_flux(p),x,y,z)
        ! Subtract each face before forming its divergence, retaining resolved
        ! changes beneath a large old bulk transfer. Identity remains exact.
        if(any(change(2:5)/=0)) work%fields(p)%values(2:5,x,y,z)= &
          history%forcing%fields(p)%values(2:5,x,y,z)+change(2:5)/volume/history%interval
        work%fields(p)%values(1,x,y,z)=history%mass_source%fields(p)%values(1,x,y,z)+ &
          new_div(1)/volume/history%interval
        work%fields(p)%values(6:,x,y,z)=history%mass_source%fields(p)%values(6:,x,y,z)+ &
          new_div(6:)/volume/history%interval
        held_change=held_mol_face_change_3d(g,scalar_flux(p),history%internal_energy_flux(p),x,y,z)
        if(held_change(1)/=0) held_work%fields(p)%values(1,x,y,z)= &
          history%internal_energy_spatial_rate%fields(p)%values(1,x,y,z)+held_change(1)/volume/history%interval
      enddo
      enddo
      enddo
      if(any(.not.ieee_is_finite(work%fields(p)%values)).or. &
        any(.not.ieee_is_finite(held_work%fields(p)%values))) return
      end associate
    enddo
    rates=work;scalar_rates=held_work;ok=.true.
  end subroutine reconstruct_reactive_eb_mol_held_interval_3d

  function held_mol_face_impulse_3d(g,f,x,y,z) result(value)
    type(eb_geometry_3d), intent(in) :: g
    type(eb_integrated_flux_3d), intent(in) :: f
    integer, intent(in) :: x,y,z
    real(dp) :: value(size(f%x,1))
    value=g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*f%x(:,x-1,y,z)-g%x_face_fraction(x,y,z)*f%x(:,x,y,z))+ &
      g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*f%y(:,x,y-1,z)-g%y_face_fraction(x,y,z)*f%y(:,x,y,z))+ &
      g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*f%z(:,x,y,z-1)-g%z_face_fraction(x,y,z)*f%z(:,x,y,z))
  end function held_mol_face_impulse_3d

  function held_mol_face_change_3d(g,f,old,x,y,z) result(value)
    type(eb_geometry_3d), intent(in) :: g
    type(eb_integrated_flux_3d), intent(in) :: f,old
    integer, intent(in) :: x,y,z
    real(dp) :: value(size(f%x,1))
    value=g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*(f%x(:,x-1,y,z)-old%x(:,x-1,y,z))- &
      g%x_face_fraction(x,y,z)*(f%x(:,x,y,z)-old%x(:,x,y,z)))+ &
      g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*(f%y(:,x,y-1,z)-old%y(:,x,y-1,z))- &
      g%y_face_fraction(x,y,z)*(f%y(:,x,y,z)-old%y(:,x,y,z)))+ &
      g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*(f%z(:,x,y,z-1)-old%z(:,x,y,z-1))- &
      g%z_face_fraction(x,y,z)*(f%z(:,x,y,z)-old%z(:,x,y,z)))
  end function held_mol_face_change_3d

  ! Replace only open-face impulses. Retain the native non-face spatial
  ! terms (including pressure work) without subtracting two bulk rates.
  ! Gas and the independent held scalar must use the same spatial family.
  subroutine blend_native_mol_face_pair_3d(g,dt,weight,low,held_low,rate,faces,held_rate,held_faces,good, &
      mass_species_from_faces)
    type(eb_geometry_3d), intent(in) :: g
    real(dp), intent(in) :: dt,weight
    type(eb_integrated_flux_3d), intent(in) :: low,held_low
    real(dp), intent(inout) :: rate(:,:,:,:),held_rate(:,:,:,:)
    type(eb_integrated_flux_3d), intent(inout) :: faces,held_faces
    logical, intent(out) :: good
    logical, optional, intent(in) :: mass_species_from_faces
    type(eb_integrated_flux_3d) :: blended,held_blended
    real(dp) :: volume,scalar_change(1)
    real(dp) :: gas_impulse(size(rate,1))
    integer :: x,y,z,nv
    good=.false.;nv=size(rate,1)
    if(.not.ieee_is_finite(dt).or..not.ieee_is_finite(weight)) return
    if(dt<=0.or.weight<0.or.weight>1) return
    if(any(shape(rate)/=[nv,g%nx,g%ny,g%nz]).or. &
      any(shape(held_rate)/=[1,g%nx,g%ny,g%nz])) return
    if(.not.valid_eb_integrated_flux_3d(faces,g,nv).or. &
      .not.valid_eb_integrated_flux_3d(held_faces,g,1)) return
    if(weight==1) then
      good=.true.;return
    endif
    if(.not.valid_eb_integrated_flux_3d(low,g,nv).or. &
      .not.valid_eb_integrated_flux_3d(held_low,g,1)) return
    blended=low;held_blended=held_low
    if(weight/=0) then
      blended%x=(1-weight)*low%x+weight*faces%x
      blended%y=(1-weight)*low%y+weight*faces%y
      blended%z=(1-weight)*low%z+weight*faces%z
      held_blended%x=(1-weight)*held_low%x+weight*held_faces%x
      held_blended%y=(1-weight)*held_low%y+weight*held_faces%y
      held_blended%z=(1-weight)*held_low%z+weight*held_faces%z
    endif
    do z=1,g%nz
    do y=1,g%ny
    do x=1,g%nx
      if(g%volume_fraction(x,y,z)<=0) cycle
      volume=g%dx*g%dy*g%dz*g%volume_fraction(x,y,z)
      rate(:,x,y,z)=rate(:,x,y,z)+held_mol_face_change_3d(g,blended,faces,x,y,z)/volume/dt
      if(present(mass_species_from_faces)) then
        if(mass_species_from_faces) then
          ! Built-in spatial mass/species have no volume source. Assemble
          ! their accepted divergence directly, rather than cancelling the
          ! old RHS with an opposing face delta and creating absent traces.
          gas_impulse=g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*blended%x(:,x-1,y,z)- &
            g%x_face_fraction(x,y,z)*blended%x(:,x,y,z))+ &
            g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*blended%y(:,x,y-1,z)- &
            g%y_face_fraction(x,y,z)*blended%y(:,x,y,z))+ &
            g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*blended%z(:,x,y,z-1)- &
            g%z_face_fraction(x,y,z)*blended%z(:,x,y,z))
          rate(1,x,y,z)=gas_impulse(1)/volume/dt
          rate(6:,x,y,z)=gas_impulse(6:)/volume/dt
        endif
      endif
      scalar_change=held_mol_face_change_3d(g,held_blended,held_faces,x,y,z)
      held_rate(1,x,y,z)=held_rate(1,x,y,z)+scalar_change(1)/volume/dt
    enddo
    enddo
    enddo
    faces=blended;held_faces=held_blended
    call post_eb_spatial_density_rate_3d(g,rate,good)
    if(good) good=all(ieee_is_finite(held_rate))
  end subroutine blend_native_mol_face_pair_3d

  logical function valid_native_mol_clock_profile_3d(starts,periods,depth,time,interval,levels) result(good)
    real(dp), allocatable, intent(in) :: starts(:),periods(:)
    real(dp), intent(in) :: time,interval
    integer, intent(in) :: depth,levels
    good=.false.
    if(depth<0.or.depth>=levels) return
    if(.not.all(ieee_is_finite([time,interval]))) return
    if(.not.allocated(starts).or..not.allocated(periods)) return
    if(size(starts)/=levels.or.size(periods)/=levels.or.lbound(starts,1)/=1.or.lbound(periods,1)/=1) return
    if(.not.all(ieee_is_finite(starts)).or..not.all(ieee_is_finite(periods))) return
    if(any(starts<0).or.any(periods<=0)) return
    if(starts(depth+1)/=time.or.periods(depth+1)/=interval) return
    good=.true.
  end function valid_native_mol_clock_profile_3d

  logical function valid_held_mol_record_3d(tree,history,nv,frame) result(good)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(eb_mol_held_interval_3d), intent(in) :: history
    integer, intent(in) :: nv
    type(rotating_frame_3d), optional, intent(in) :: frame
    integer(int64), allocatable :: origin(:,:),domain(:,:),initial_origin(:,:),initial_domain(:,:)
    integer, allocatable :: ratios(:),initial_ratios(:)
    type(eb_geometry_3d), allocatable :: funding_geometry(:)
    integer :: p,x,y,z,n(3),stage,allocation_status
    logical :: valid
    good=.false.
    if(history%auxiliary_has_internal_energy.and..not.history%has_auxiliary_source) return
    if(history%routed_auxiliary) then
      if(.not.history%has_auxiliary_source) return
      if(.not.valid_held_mol_field_3d(tree,history%routed_source, &
        eb_auxiliary_nvar(nv,history%auxiliary_has_internal_energy))) return
    endif
    if(.not.valid_tree(tree,nv).or..not.valid_tree(history%initial,nv)) return
    if(history%depth<0.or.history%depth>=tree%topology%level_count()) return
    if(allocated(history%clock_source_funding_context%pending).neqv. &
      allocated(history%clock_source_funding_context%open)) return
    if(allocated(history%clock_source_funding_context%pending)) then
      if(.not.history%has_source_funding) return
      allocate(funding_geometry(size(tree%nodes)),stat=allocation_status)
      if(allocation_status/=0) return
      do p=1,size(tree%nodes)
        funding_geometry(p)=tree%nodes(p)%geometry
      enddo
      if(.not.valid_eb_held_funding_context_3d(tree%topology,funding_geometry,nv, &
        [( .true.,p=1,size(tree%nodes))],history%clock_source_funding_context,history%depth)) return
    else
      if(allocated(history%clock_source_funding_context%topology%nodes).or. &
        allocated(history%clock_source_funding_context%geometry)) return
    endif
    if(.not.all(ieee_is_finite([history%time,history%interval,history%forcing_energy_correction, &
      history%numerical_energy_correction,history%clock_energy_correction,history%internal_energy_projection_change, &
      history%internal_energy_closure_change,history%initial_internal_energy_change, &
      history%held_clock_start_change,history%held_clock_end_change, &
      history%held_clock_projection_before,history%held_clock_reset_before, &
      history%held_clock_projection_after,history%held_clock_reset_after]))) return
    if(history%held_clock_prefixes_valid.and..not.history%held_clock_events_closed) return
    if(history%time<0.or.history%interval<=0) return
    if(history%native_iterations<0.or.any(.not.ieee_is_finite(history%native_force))) return
    if(history%native_iterations==0) then
      if(allocated(history%external_stages)) return
      if(allocated(history%auxiliary_stages).or.allocated(history%provider_stages).or.any(history%native_providers)) return
    else
      if(history%prepared_redistribution) return
      if(.not.allocated(history%external_stages)) return
      if(lbound(history%external_stages,1)/=0.or.ubound(history%external_stages,1)/=history%native_iterations) return
      do stage=0,history%native_iterations
        if(.not.valid_held_mol_field_3d(tree,history%external_stages(stage),nv)) return
      enddo
      if(history%has_auxiliary_source.neqv.allocated(history%auxiliary_stages)) return
      if(allocated(history%auxiliary_stages)) then
        if(lbound(history%auxiliary_stages,1)/=0.or.ubound(history%auxiliary_stages,1)/=history%native_iterations) return
        do stage=0,history%native_iterations
          if(.not.valid_held_mol_field_3d(tree,history%auxiliary_stages(stage), &
            eb_auxiliary_nvar(nv,history%auxiliary_has_internal_energy))) return
        enddo
      endif
    endif
    if(history%has_source_funding) then
      if(history%prepared_redistribution) return
      if(.not.valid_held_mol_field_3d(tree,history%source_funding_stock,nv+1).or. &
        .not.valid_held_mol_field_3d(tree,history%source_funding_impulse,nv+1)) return
      do p=1,size(tree%nodes)
        if(tree%topology%nodes(p)%level==history%depth) then
          if(any(history%source_funding_stock%fields(p)%values/=0).or. &
            any(history%source_funding_impulse%fields(p)%values/=0)) return
        endif
        n=tree%topology%nodes(p)%shape
        do z=1,n(3)
        do y=1,n(2)
        do x=1,n(1)
          if(.not.covered(tree%topology,p,[x,y,z]).and.tree%nodes(p)%geometry%volume_fraction(x,y,z)>0) cycle
          if(any(history%source_funding_stock%fields(p)%values(:,x,y,z)/=0).or. &
            any(history%source_funding_impulse%fields(p)%values(:,x,y,z)/=0)) return
        enddo
        enddo
        enddo
      enddo
    else
      if(allocated(history%source_funding_stock%fields).or.allocated(history%source_funding_impulse%fields)) return
    endif
    if(history%has_frame.neqv.present(frame)) return
    if(present(frame)) then
      if(.not.frame%is_valid().or..not.history%frame%is_valid()) return
      if(frame%axis/=history%frame%axis.or.frame%angular_speed/=history%frame%angular_speed.or. &
        any(frame%axis_point/=history%frame%axis_point)) return
    endif
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if(.not.valid) return
    call amr_tree_integer_layout_3d(history%initial%topology,initial_origin,initial_domain,initial_ratios,valid)
    if(.not.valid) return
    if(any(shape(origin)/=shape(initial_origin)).or.any(shape(domain)/=shape(initial_domain))) return
    if(any(origin/=initial_origin).or.any(domain/=initial_domain).or.any(ratios/=initial_ratios)) return
    if(.not.valid_held_mol_field_3d(tree,history%forcing,nv).or. &
      .not.valid_held_mol_field_3d(tree,history%mass_source,nv).or. &
      .not.valid_held_mol_field_3d(tree,history%initial_internal_energy,1).or. &
      .not.valid_held_mol_field_3d(tree,history%internal_energy_spatial_rate,1).or. &
      .not.valid_held_mol_field_3d(tree,history%initial_reaction_source,nv-4)) return
    if(.not.allocated(history%spatial_flux).or..not.allocated(history%internal_energy_flux)) return
    if(size(history%spatial_flux)/=size(tree%nodes).or.size(history%internal_energy_flux)/=size(tree%nodes)) return
    do p=1,size(tree%nodes)
      if(tree%topology%nodes(p)%parent/=history%initial%topology%nodes(p)%parent) return
      if(any([tree%nodes(p)%geometry%dx,tree%nodes(p)%geometry%dy,tree%nodes(p)%geometry%dz, &
        tree%nodes(p)%geometry%x_lower,tree%nodes(p)%geometry%y_lower,tree%nodes(p)%geometry%z_lower]/= &
        [history%initial%nodes(p)%geometry%dx,history%initial%nodes(p)%geometry%dy,history%initial%nodes(p)%geometry%dz, &
        history%initial%nodes(p)%geometry%x_lower,history%initial%nodes(p)%geometry%y_lower, &
        history%initial%nodes(p)%geometry%z_lower])) return
      if(.not.matching_eb_metrics_3d(tree%nodes(p)%geometry,history%initial%nodes(p)%geometry)) return
      if(any(.not.ieee_is_finite(history%initial%nodes(p)%state)).or. &
        any(.not.ieee_is_finite(history%initial%nodes(p)%temperature))) return
      if(tree%topology%nodes(p)%level/=history%depth.and..not.history%has_source_funding) cycle
      if(.not.valid_eb_integrated_flux_3d(history%spatial_flux(p),tree%nodes(p)%geometry,nv).or. &
        .not.valid_eb_integrated_flux_3d(history%internal_energy_flux(p),tree%nodes(p)%geometry,1)) return
    enddo
    good=.true.
  end function valid_held_mol_record_3d

  logical function valid_mol_native_provider_cache_3d(tree,history,nv) result(good)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(eb_mol_held_interval_3d), intent(in) :: history
    integer, intent(in) :: nv
    integer :: stage,kind,p,d(3)
    logical :: expected
    good=.false.
    expected=history%spatial_diffusion_is_source.and.any(history%native_providers)
    if(allocated(history%provider_stages).neqv.expected) return
    if(.not.expected) then
      good=.true.;return
    endif
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
        if(.not.valid_held_mol_field_3d(tree,saved%provider_rate,nv)) return
        if(.not.allocated(saved%provider_faces)) return
        if(size(saved%provider_faces)/=size(tree%nodes).or.lbound(saved%provider_faces,1)/=1) return
        if(kind==1) then
          if(.not.valid_held_prepared_halo_3d(tree,saved%provider_halos,nv,history%depth)) return
        else
          if(allocated(saved%provider_halos%fields)) return
        endif
        do p=1,size(tree%nodes)
          if(tree%topology%nodes(p)%level/=history%depth) cycle
          if(.not.valid_eb_integrated_flux_3d(saved%provider_faces(p),tree%nodes(p)%geometry,nv)) return
          if(kind/=1) cycle
          d=tree%topology%nodes(p)%shape
          if(any(lbound(saved%provider_halos%fields(p)%values)/=1)) return
          if(any(saved%provider_halos%fields(p)%values(:,7:d(1)+6,7:d(2)+6,7:d(3)+6)/= &
            saved%provider_rate%fields(p)%values)) return
        enddo
        end associate
      enddo
    enddo
    good=.true.
  end function valid_mol_native_provider_cache_3d

  logical function valid_held_mol_field_3d(tree,field,nv) result(good)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: field
    integer, intent(in) :: nv
    integer :: p
    good=.false.
    if(.not.field%is_valid(tree%topology)) return
    do p=1,size(tree%nodes)
      if(size(field%fields(p)%values,1)/=nv.or.any(.not.ieee_is_finite(field%fields(p)%values))) return
    enddo
    good=.true.
  end function valid_held_mol_field_3d

  ! Reconstruct the fixed-reference *registration* from admitted source
  ! stages. The accepted old/final state-update registers are a separate
  ! family; this routine neither changes that state nor calls chemistry.
  subroutine reconstruct_reactive_eb_mol_held_registration_3d(tree,history,registers,held_registers,ok,frame)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(eb_mol_held_interval_3d), intent(in) :: history
    type(state_redist_register_3d), allocatable, intent(inout) :: registers(:),held_registers(:)
    logical, intent(out) :: ok
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(state_redist_register_3d), allocatable :: pending(:),held_pending(:)
    integer :: nv,steps,stage,p,d(3),ios
    logical :: good
    ok=.false.
    if(.not.allocated(tree%nodes)) return
    if(size(tree%nodes)<1) return
    if(.not.allocated(tree%nodes(1)%state)) return
    nv=size(tree%nodes(1)%state,1)
    if(.not.allocated(history%prepared_stages)) return
    steps=ubound(history%prepared_stages,1)
    if(steps<1.or.lbound(history%prepared_stages,1)/=0) return
    if(.not.valid_held_prepared_mol_record_3d(tree,history,nv,history%depth,history%interval,steps, &
      history%prepared_stages(0)%target_volume_fraction,history%prepared_stages(0)%max_order,history%policy, &
      frame,history%has_auxiliary_source)) return
    allocate(pending(size(tree%nodes)),held_pending(size(tree%nodes)),stat=ios)
    if(ios/=0) return
    do stage=0,steps
      associate(saved=>history%prepared_stages(stage))
      if(.not.allocated(saved%registers).or..not.allocated(saved%scalar_registers)) return
      if(size(saved%registers)/=size(tree%nodes).or.size(saved%scalar_registers)/=size(tree%nodes)) return
      if(lbound(saved%registers,1)/=1.or.lbound(saved%scalar_registers,1)/=1) return
      do p=1,size(tree%nodes)
        if(tree%topology%nodes(p)%level/=history%depth) cycle
        d=tree%topology%nodes(p)%shape+12
        call accumulate_state_redist_mol_registration_stage_3d(pending(p),saved%registers(p),stage,steps,nv,d,good)
        if(.not.good) return
        call accumulate_state_redist_mol_registration_stage_3d(held_pending(p),saved%scalar_registers(p), &
          stage,steps,1,d,good)
        if(.not.good) return
      enddo
      end associate
    enddo
    call move_alloc(pending,registers);call move_alloc(held_pending,held_registers)
    ok=.true.
  end subroutine reconstruct_reactive_eb_mol_held_registration_3d

  logical function valid_held_prepared_mol_record_3d(tree,history,nv,depth,interval,iterations, &
      redist_target,redist_order,policy,frame,auxiliary_source_supplied) result(good)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(eb_mol_held_interval_3d), intent(in) :: history
    integer, intent(in) :: nv,depth,iterations,redist_order
    real(dp), intent(in) :: interval,redist_target
    type(dual_energy_policy), intent(in) :: policy
    type(rotating_frame_3d), optional, intent(in) :: frame
    integer :: stage
    logical, optional, intent(in) :: auxiliary_source_supplied
    logical :: has_auxiliary
    good=.false.;has_auxiliary=.false.
    if(present(auxiliary_source_supplied)) has_auxiliary=auxiliary_source_supplied
    if(history%separate_momentum) return
    if(.not.history%prepared_redistribution.or..not.history%accepted_coarse_stages) return
    if(history%auxiliary_has_internal_energy.and..not.history%has_auxiliary_source) return
    if(history%has_auxiliary_source.neqv.has_auxiliary) return
    if(history%depth/=depth.or.history%interval/=interval) return
    if(.not.valid_held_mol_record_3d(tree,history,nv,frame)) return
    if(history%policy%eta2/=policy%eta2.or. &
      (history%policy%allow_small_energy.neqv.policy%allow_small_energy).or. &
      (history%policy%allow_negative_energy.neqv.policy%allow_negative_energy).or. &
      (history%policy%update_total_from_internal.neqv.policy%update_total_from_internal)) return
    if(.not.allocated(history%prepared_stages)) return
    if(lbound(history%prepared_stages,1)/=0.or.ubound(history%prepared_stages,1)/=iterations) return
    do stage=0,iterations
      associate(saved=>history%prepared_stages(stage))
      if(saved%source_stage/=stage.or.saved%depth/=depth.or.saved%interval/=interval.or. &
        saved%max_order/=redist_order.or.saved%target_volume_fraction/=redist_target.or.saved%register_scale/=.5_dp) return
      if(.not.ieee_is_finite(saved%spatial_weight).or.saved%spatial_weight<0.or.saved%spatial_weight>1) return
      if(.not.valid_held_mol_field_3d(tree,saved%external_rate,nv).or. &
        .not.valid_held_mol_field_3d(tree,saved%auxiliary_rate,eb_auxiliary_nvar(nv, &
      history%auxiliary_has_internal_energy))) return
      if(.not.valid_held_prepared_halo_3d(tree,saved%state_halos,nv,depth).or. &
        .not.valid_held_prepared_halo_3d(tree,saved%raw,nv,depth).or. &
        .not.valid_held_prepared_halo_3d(tree,saved%scalar_state_halos,1,depth).or. &
        .not.valid_held_prepared_halo_3d(tree,saved%scalar_raw,1,depth)) return
      end associate
    enddo
    good=.true.
  end function valid_held_prepared_mol_record_3d

  logical function valid_held_prepared_halo_3d(tree,field,components,depth) result(good)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: field
    integer, intent(in) :: components,depth
    integer :: p,d(3)
    good=.false.
    if(.not.allocated(field%fields)) return
    if(size(field%fields)/=size(tree%nodes)) return
    do p=1,size(tree%nodes)
      if(tree%topology%nodes(p)%level/=depth) cycle
      if(.not.allocated(field%fields(p)%values)) return
      d=tree%topology%nodes(p)%shape
      if(any(shape(field%fields(p)%values)/=[components,d+12])) return
      if(any(.not.ieee_is_finite(field%fields(p)%values))) return
    enddo
    good=.true.
  end function valid_held_prepared_halo_3d

  ! Native ordinary held replay. The scalar reconstruction is archived, while
  ! forced chemistry receives projected gas forcing and the actual old UEINT.
  ! Every field, heat payload, reset charge and replacement record is private
  ! until all selected cells pass the unchanged reaction/EOS gates.
  subroutine replay_reactive_eb_mol_held_interval_3d(species,reactions,tree,stored,reaction_source,prefix,history, &
      flux,scalar_flux,chemical_change,rtol,atol,ok,evaluate_reaction,frame,heat_release,provider_refused, &
      source_funding_stock,source_funding_impulse)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(eb_reactive_tree_3d), intent(inout) :: tree
    type(amr_tree_state_3d), intent(inout) :: stored,reaction_source
    real(dp), intent(inout) :: prefix,chemical_change(:)
    type(eb_mol_held_interval_3d), intent(inout) :: history
    type(eb_integrated_flux_3d), intent(in) :: flux(:),scalar_flux(:)
    real(dp), intent(in) :: rtol,atol
    logical, intent(out) :: ok
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d), optional, intent(inout) :: heat_release,source_funding_impulse
    type(amr_tree_state_3d), optional, intent(in) :: source_funding_stock
    logical, optional, intent(out) :: provider_refused
    type(eb_reactive_tree_3d) :: old,candidate,ledger
    type(amr_tree_state_3d) :: held,old_held,ir,rates,scalar_rates,u,t,reset_delta
    type(amr_tree_state_3d), allocatable :: heat
    type(eb_mol_held_interval_3d) :: record
    procedure(forced_reactor_integrator), pointer :: integrator
    real(dp) :: initial_reset,final_reset,correction,amount(reactive_nvar(size(species)))
    integer :: nv,p,x,y,z,n(3)
    real(dp) :: projection(1),transported
    real(dp) :: initial_scalar_reset
    type(amr_tree_state_3d) :: projection_delta
    logical :: good,refused
    ok=.false.;nv=reactive_nvar(size(species));refused=.false.;nullify(integrator)
    if(present(provider_refused)) provider_refused=.false.
    if(.not.ieee_is_finite(prefix).or.size(chemical_change)/=nv) return
    if(.not.valid_held_mol_field_3d(tree,stored,1).or. &
      .not.valid_held_mol_field_3d(tree,reaction_source,nv-4)) return
    if(present(heat_release)) then
      if(.not.valid_held_mol_field_3d(tree,heat_release,1)) return
      allocate(heat,source=heat_release)
    endif
    call reconstruct_reactive_eb_mol_held_interval_3d(species,tree,history,flux,scalar_flux,rates,scalar_rates,good,frame)
    if(.not.good) return
    if(present(source_funding_stock).neqv.present(source_funding_impulse)) return
    if(present(source_funding_stock).neqv.history%has_source_funding) return
    if(present(source_funding_stock)) then
      if(.not.valid_held_mol_field_3d(tree,source_funding_stock,nv+1)) return
      do p=1,size(tree%nodes)
        if(any(source_funding_stock%fields(p)%values/=history%source_funding_stock%fields(p)%values)) return
        if(tree%topology%nodes(p)%level==history%depth) cycle
        if(any(tree%nodes(p)%state/=history%initial%nodes(p)%state).or. &
          any(stored%fields(p)%values/=history%initial_internal_energy%fields(p)%values)) return
        if(any(flux(p)%x/=history%spatial_flux(p)%x).or.any(flux(p)%y/=history%spatial_flux(p)%y).or. &
          any(flux(p)%z/=history%spatial_flux(p)%z).or. &
          any(scalar_flux(p)%x/=history%internal_energy_flux(p)%x).or. &
          any(scalar_flux(p)%y/=history%internal_energy_flux(p)%y).or. &
          any(scalar_flux(p)%z/=history%internal_energy_flux(p)%z)) return
      enddo
    endif
    old=tree;old_held=stored;ir=reaction_source;record=history
    call reset_held_eb_tree_family_3d(species,old,old_held,reset_delta,history%policy,history%depth, &
      initial_reset,good,frame, &
      internal_energy_change=initial_scalar_reset,reset_level=-2)
    if(.not.good) return
    record%clock_initial=old;record%clock_internal_energy=old_held;record%clock_reaction_source=ir
    record%clock_energy_correction=prefix+initial_reset
    if(present(heat_release)) then
      record%clock_heat=heat
    else
      if(allocated(record%clock_heat%fields)) deallocate(record%clock_heat%fields)
    endif
    record%initial_reaction_source=ir
    u=rates;t=old_held;held=old_held
    do p=1,size(tree%nodes)
      u%fields(p)%values=old%nodes(p)%state
      t%fields(p)%values(1,:,:,:)=old%nodes(p)%temperature
    enddo
    if(present(evaluate_reaction)) integrator=>classified_reactor
    if(associated(integrator)) then
      call react_sdc_eb_tree_3d(species,reactions,old,rates,history%interval,rtol,atol,u,t,ir,good, &
        target_depth=history%depth,evaluate_reaction=integrator,frame=frame,heat_release=heat, &
        internal_energy=held,old_internal_energy=old_held)
    else
      call react_sdc_eb_tree_3d(species,reactions,old,rates,history%interval,rtol,atol,u,t,ir,good, &
        target_depth=history%depth,frame=frame,heat_release=heat,internal_energy=held,old_internal_energy=old_held)
    endif
    if(present(provider_refused)) provider_refused=refused
    if(.not.good.or.refused) return
    candidate=old;ledger=old
    do p=1,size(tree%nodes)
      candidate%nodes(p)%state=u%fields(p)%values;candidate%nodes(p)%temperature=t%fields(p)%values(1,:,:,:)
    enddo
    call reset_held_eb_tree_family_3d(species,candidate,held,reset_delta,history%policy,history%depth, &
      final_reset,good,frame)
    if(.not.good) return
    do p=1,size(tree%nodes)
      ledger%nodes(p)%state=0
      if(tree%topology%nodes(p)%level/=history%depth) cycle
      ledger%nodes(p)%state=candidate%nodes(p)%state-(old%nodes(p)%state+history%interval*rates%fields(p)%values)
    enddo
    call integrate_reactive_eb_tree_3d(ledger,amount,good)
    if(.not.good) return
    amount(5)=amount(5)-final_reset
    correction=prefix+initial_reset+history%forcing_energy_correction+final_reset
    if(.not.ieee_is_finite(correction).or.any(.not.ieee_is_finite(amount))) return
    record%initial=old;record%initial_internal_energy=old_held;record%forcing=rates
    record%local_tree=candidate;record%local_internal_energy=held;record%reaction_source=ir
    record%internal_energy_spatial_rate=scalar_rates;record%chemical_change=amount
    record%spatial_flux=flux;record%internal_energy_flux=scalar_flux
    record%numerical_energy_correction=initial_reset+history%forcing_energy_correction+final_reset
    projection_delta=old_held
    do p=1,size(tree%nodes)
      projection_delta%fields(p)%values=0
      if(tree%topology%nodes(p)%level/=history%depth) cycle
      n=tree%topology%nodes(p)%shape
      do z=1,n(3)
      do y=1,n(2)
      do x=1,n(1)
        if(old%nodes(p)%geometry%volume_fraction(x,y,z)<=0.or.covered(tree%topology,p,[x,y,z])) cycle
        transported=old_held%fields(p)%values(1,x,y,z)+history%interval*scalar_rates%fields(p)%values(1,x,y,z)
        call compute_reactive_eb_held_projection(stored%fields(p)%values(1,x,y,z), &
          old_held%fields(p)%values(1,x,y,z),transported,held%fields(p)%values(1,x,y,z), &
          projection_delta%fields(p)%values(1,x,y,z),good)
        if(.not.good) return
      enddo
      enddo
      enddo
      ledger%nodes(p)%state=projection_delta%fields(p)%values
    enddo
    ! Inactive nodes still need a scalar, zero-valued integral payload.
    do p=1,size(tree%nodes)
      ledger%nodes(p)%state=projection_delta%fields(p)%values
    enddo
    projection=0
    call integrate_reactive_eb_tree_3d(ledger,projection,good)
    if(.not.good) return
    record%internal_energy_projection_change=projection(1)
    record%initial_internal_energy_change=initial_scalar_reset
    record%held_clock_start_change=initial_scalar_reset;record%held_clock_end_change=0
    ! This replacement ends at its local solve; any following hierarchy
    ! reset is owned by the surrounding clock group, not the previous record.
    record%internal_energy_closure_change=0
    record%held_clock_events_closed=.false.;record%held_clock_prefixes_valid=.false.
    record%held_clock_projection_before=0;record%held_clock_reset_before=0
    record%held_clock_projection_after=0;record%held_clock_reset_after=0
    record%replayed=.true.
    if(present(source_funding_impulse)) source_funding_impulse=record%source_funding_impulse
    tree=candidate;stored=held;reaction_source=ir;prefix=correction;chemical_change=amount;history=record
    if(present(heat_release)) heat_release=heat
    ok=.true.
  contains
    subroutine classified_reactor(species,reactions,energy,source_species,source_energy,duration,rtol,atol, &
        density,temperature,steps,ok,maximum_steps)
      type(nasa7_species), intent(in) :: species(:)
      type(elementary_reaction), intent(in) :: reactions(:)
      real(dp), intent(in) :: energy,source_species(:),source_energy,duration,rtol,atol
      real(dp), intent(inout) :: density(:),temperature
      integer, intent(out) :: steps
      logical, intent(out) :: ok
      integer, optional, intent(in) :: maximum_steps
      real(dp) :: probe(size(density)),probe_t
      integer :: probe_steps
      logical :: native_ok
      probe=density;probe_t=temperature
      call evaluate_reaction(species,reactions,energy,source_species,source_energy,duration,rtol,atol, &
        density,temperature,steps,ok,maximum_steps)
      if(ok) return
      call advance_forced_reactor_interval(species,reactions,energy,source_species,source_energy,duration,rtol,atol, &
        probe,probe_t,probe_steps,native_ok,maximum_steps=maximum_steps)
      refused=refused.or.native_ok
    end subroutine classified_reactor
  end subroutine replay_reactive_eb_mol_held_interval_3d

  subroutine allocate_held_mol_faces_3d(tree,nv,faces)
    type(eb_reactive_tree_3d), intent(in) :: tree
    integer, intent(in) :: nv
    type(eb_integrated_flux_3d), allocatable, intent(out) :: faces(:)
    integer :: p,n(3)
    allocate(faces(size(tree%nodes)))
    do p=1,size(tree%nodes)
      n=tree%topology%nodes(p)%shape
      allocate(faces(p)%x(nv,0:n(1),n(2),n(3)),faces(p)%y(nv,n(1),0:n(2),n(3)), &
        faces(p)%z(nv,n(1),n(2),0:n(3)),source=0._dp)
    enddo
  end subroutine allocate_held_mol_faces_3d

  ! Restore and integrate a complete physical parent window, rather than repair
  ! a chemically consumed coarse endpoint with a donor-only state increment.
  ! All gas/UEINT faces use the same weight and clock durations. The low anchor
  ! is the first old spatial family on each level in this window, avoiding a
  ! later fine old family already contaminated by the discarded parent solve.
  subroutine replay_reactive_eb_mol_held_group_3d(species,reactions,tree,stored,reaction_source,prefix,history, &
      weight,flux,scalar_flux,chemical_change,rtol,atol,ok,boundary_faces,evaluate_reaction,frame,heat_release, &
      provider_refused)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(eb_reactive_tree_3d), intent(inout) :: tree
    type(amr_tree_state_3d), intent(inout) :: stored,reaction_source
    real(dp), intent(inout) :: prefix,chemical_change(:)
    type(eb_mol_held_interval_3d), intent(inout) :: history(:)
    real(dp), intent(in) :: weight,rtol,atol
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:),scalar_flux(:)
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d), optional, intent(inout) :: heat_release
    logical, optional, intent(out) :: provider_refused
    type :: paired_family
      type(eb_integrated_flux_3d), allocatable :: faces(:),scalar_faces(:)
    end type paired_family
    type(paired_family), allocatable :: family(:)
    type(eb_mol_held_interval_3d), allocatable :: records(:)
    type(eb_reactive_tree_3d) :: work
    type(amr_tree_state_3d) :: held,ir,increment,scratch,reset_delta,packed_gas
    type(amr_tree_state_3d), allocatable :: routed_pending(:)
    type(amr_tree_state_3d), allocatable :: heat
    type(eb_integrated_flux_3d), allocatable :: raw(:),canonical(:),scalar_raw(:),result(:),scalar_result(:)
    real(dp), allocatable :: cursor(:)
    integer, allocatable :: reference(:)
    real(dp) :: scale,tolerance,start,ending,correction,closure_reset,amount(reactive_nvar(size(species)))
    real(dp) :: chemical(reactive_nvar(size(species))),scalar_reset,routed_scalar_reset,accepted_start
    integer :: closing_record
    integer :: nv,levels,depth,r,q,p,lev,next_depth,closure,bucket,target
    logical :: good,refused
    ok=.false.;nv=reactive_nvar(size(species))
    if(present(provider_refused)) provider_refused=.false.
    if(size(history)<1.or.size(chemical_change)/=nv.or..not.ieee_is_finite(prefix)) return
    if(.not.ieee_is_finite(weight).or.weight<0.or.weight>1) return
    levels=tree%topology%level_count();depth=history(1)%depth
    if(depth<0.or.depth>=levels) return
    allocate(reference(0:levels-1),cursor(0:levels-1),family(size(history)))
    reference=0;start=history(1)%time;ending=start+history(1)%interval;cursor=start
    if(.not.ieee_is_finite(ending).or.ending<=start) return
    tolerance=32*epsilon(1._dp)*max(abs(ending),tiny(1._dp))
    do r=1,size(history)
      if(.not.valid_held_mol_record_3d(tree,history(r),nv,frame)) return
      if(history(r)%auxiliary_has_internal_energy.neqv.history(1)%auxiliary_has_internal_energy) return
      if(history(r)%routed_auxiliary.neqv.history(1)%routed_auxiliary) return
      if(history(r)%policy%eta2/=history(1)%policy%eta2.or. &
        (history(r)%policy%allow_small_energy.neqv.history(1)%policy%allow_small_energy).or. &
        (history(r)%policy%allow_negative_energy.neqv.history(1)%policy%allow_negative_energy).or. &
        (history(r)%policy%update_total_from_internal.neqv.history(1)%policy%update_total_from_internal)) return
      lev=history(r)%depth
      if(lev<depth) return
      if(abs(history(r)%time-cursor(lev))>tolerance) return
      cursor(lev)=history(r)%time+history(r)%interval
      if(.not.ieee_is_finite(cursor(lev)).or.cursor(lev)>ending+tolerance) return
      if(reference(lev)==0) reference(lev)=r
      if(r>1) then
        if(lev>history(r-1)%depth+1) return
      endif
      if(.not.allocated(history(r)%old_flux).or..not.allocated(history(r)%old_internal_energy_flux)) return
      if(size(history(r)%old_flux)/=size(tree%nodes).or. &
        size(history(r)%old_internal_energy_flux)/=size(tree%nodes)) return
      do p=1,size(tree%nodes)
        if(tree%topology%nodes(p)%level/=lev) cycle
        if(.not.valid_eb_integrated_flux_3d(history(r)%old_flux(p),tree%nodes(p)%geometry,nv).or. &
          .not.valid_eb_integrated_flux_3d(history(r)%old_internal_energy_flux(p),tree%nodes(p)%geometry,1)) return
      enddo
    enddo
    if(any(reference(depth:)<1).or.any(abs(cursor(depth:)-ending)>tolerance)) return
    if(.not.valid_tree(history(1)%clock_initial,nv)) return
    ! Construct all canonical face families before admitting any reaction.
    do r=size(history),1,-1
      lev=history(r)%depth;q=reference(lev);scale=history(r)%interval/history(q)%interval
      call allocate_held_mol_faces_3d(tree,nv,raw);call allocate_held_mol_faces_3d(tree,1,scalar_raw)
      do p=1,size(tree%nodes)
        if(tree%topology%nodes(p)%level/=lev) cycle
        if(weight==0) then
          raw(p)%x=scale*history(q)%old_flux(p)%x
          raw(p)%y=scale*history(q)%old_flux(p)%y
          raw(p)%z=scale*history(q)%old_flux(p)%z
          scalar_raw(p)%x=scale*history(q)%old_internal_energy_flux(p)%x
          scalar_raw(p)%y=scale*history(q)%old_internal_energy_flux(p)%y
          scalar_raw(p)%z=scale*history(q)%old_internal_energy_flux(p)%z
        else
          raw(p)%x=(1-weight)*scale*history(q)%old_flux(p)%x+weight*history(r)%spatial_flux(p)%x
          raw(p)%y=(1-weight)*scale*history(q)%old_flux(p)%y+weight*history(r)%spatial_flux(p)%y
          raw(p)%z=(1-weight)*scale*history(q)%old_flux(p)%z+weight*history(r)%spatial_flux(p)%z
          scalar_raw(p)%x=(1-weight)*scale*history(q)%old_internal_energy_flux(p)%x+ &
            weight*history(r)%internal_energy_flux(p)%x
          scalar_raw(p)%y=(1-weight)*scale*history(q)%old_internal_energy_flux(p)%y+ &
            weight*history(r)%internal_energy_flux(p)%y
          scalar_raw(p)%z=(1-weight)*scale*history(q)%old_internal_energy_flux(p)%z+ &
            weight*history(r)%internal_energy_flux(p)%z
        endif
      enddo
      do q=r+1,size(history)
        if(history(q)%depth/=lev+1) cycle
        if(history(q)%time<history(r)%time-tolerance.or. &
          history(q)%time+history(q)%interval>history(r)%time+history(r)%interval+tolerance) cycle
        do p=1,size(tree%nodes)
          if(tree%topology%nodes(p)%level/=lev+1) cycle
          raw(p)%x=raw(p)%x+family(q)%faces(p)%x
          raw(p)%y=raw(p)%y+family(q)%faces(p)%y
          raw(p)%z=raw(p)%z+family(q)%faces(p)%z
          scalar_raw(p)%x=scalar_raw(p)%x+family(q)%scalar_faces(p)%x
          scalar_raw(p)%y=scalar_raw(p)%y+family(q)%scalar_faces(p)%y
          scalar_raw(p)%z=scalar_raw(p)%z+family(q)%scalar_faces(p)%z
        enddo
      enddo
      call assemble_reactive_eb_tree_reflux_3d(species,history(1)%clock_initial,raw,increment,canonical,good, &
        lev,min(lev+1,levels-1),boundary_faces,frame)
      if(.not.good) return
      scratch=history(1)%initial_internal_energy
      do p=1,size(tree%nodes)
        scratch%fields(p)%values=0
      enddo
      call synchronize_reactive_eb_tree_scalar_flux_3d(species,history(1)%clock_initial,scratch,scalar_raw,good, &
        lev,min(lev+1,levels-1),boundary_faces)
      if(.not.good) return
      call allocate_held_mol_faces_3d(tree,nv,family(r)%faces)
      call allocate_held_mol_faces_3d(tree,1,family(r)%scalar_faces)
      do p=1,size(tree%nodes)
        if(tree%topology%nodes(p)%level/=lev) cycle
        family(r)%faces(p)=canonical(p);family(r)%scalar_faces(p)=scalar_raw(p)
      enddo
    enddo
    ! Descendant/ancestor tails and their numerical charges are discarded.
    if(.not.valid_tree(history(1)%clock_initial,nv)) return
    if(.not.valid_held_mol_field_3d(tree,history(1)%clock_internal_energy,1).or. &
      .not.valid_held_mol_field_3d(tree,history(1)%clock_reaction_source,nv-4)) return
    work=history(1)%clock_initial;held=history(1)%clock_internal_energy
    ir=history(1)%clock_reaction_source;correction=history(1)%clock_energy_correction;records=history
    if(present(heat_release)) then
      if(.not.valid_held_mol_field_3d(tree,history(1)%clock_heat,1)) return
      allocate(heat,source=history(1)%clock_heat)
    endif
    if(history(1)%routed_auxiliary) then
      allocate(routed_pending(depth:levels-1))
      do lev=depth,levels-1
        routed_pending(lev)=history(1)%routed_source
        do p=1,size(tree%nodes)
          routed_pending(lev)%fields(p)%values=0
        enddo
      enddo
    endif
    records%held_clock_start_change=0;records%held_clock_end_change=0
    records%held_clock_events_closed=.false.;records%held_clock_prefixes_valid=.false.
    chemical=0
    do r=1,size(records)
      lev=records(r)%depth
      if(allocated(routed_pending)) then
        do p=1,size(tree%nodes)
          target=tree%topology%nodes(p)%level;bucket=min(lev,target)
          if(bucket<depth) cycle
          routed_pending(bucket)%fields(p)%values=routed_pending(bucket)%fields(p)%values+ &
            records(r)%routed_source%fields(p)%values
        enddo
      endif
      call replay_reactive_eb_mol_held_interval_3d(species,reactions,work,held,ir,correction,records(r), &
        family(r)%faces,family(r)%scalar_faces,amount,rtol,atol,good, &
        evaluate_reaction,frame,heat,refused)
      if(present(provider_refused)) provider_refused=refused
      if(.not.good.or.refused) return
      ! The first canonical anchor already includes its accepted old entry
      ! repair. Keep that historical event once; later entries are remeasured.
      accepted_start=records(r)%initial_internal_energy_change
      if(r==1) then
        call compute_reactive_eb_held_projection(0._dp,history(1)%held_clock_start_change,0._dp, &
          accepted_start,scalar_reset,good)
        if(.not.good) return
        accepted_start=scalar_reset
      endif
      records(r)%held_clock_start_change=accepted_start
      chemical=chemical+amount
      call average_down_reactive_eb_tree_field_3d(work,held,good)
      if(.not.good) return
      call average_down_reactive_eb_tree_3d(species,work,good,frame)
      if(.not.good) return
      next_depth=depth
      if(r<size(records)) next_depth=records(r+1)%depth
      do closure=lev,max(next_depth,depth+1),-1
        closing_record=0
        do q=r,1,-1
          if(records(q)%depth/=closure) cycle
          closing_record=q;exit
        enddo
        if(closing_record==0) return
        scalar_reset=0
        call reset_held_eb_tree_family_3d(species,work,held,reset_delta,records(r)%policy,closure, &
          closure_reset,good,frame,leaf_only=.true.,internal_energy_change=scalar_reset)
        if(.not.good) return
        correction=correction+closure_reset
        if(allocated(routed_pending)) then
          if(allocated(packed_gas%fields)) deallocate(packed_gas%fields)
          allocate(packed_gas%fields(size(tree%nodes)))
          do p=1,size(tree%nodes)
            packed_gas%fields(p)%values=work%nodes(p)%state
          enddo
          call apply_eb_held_auxiliary_increment_3d(packed_gas,held,routed_pending(closure),good)
          if(.not.good) return
          do p=1,size(tree%nodes)
            work%nodes(p)%state=packed_gas%fields(p)%values
            routed_pending(closure)%fields(p)%values=0
          enddo
          call reset_held_eb_tree_family_3d(species,work,held,reset_delta,records(r)%policy,-1, &
            closure_reset,good,frame,leaf_only=.true.,internal_energy_change=routed_scalar_reset)
          if(.not.good) return
          call compute_reactive_eb_held_projection(0._dp,scalar_reset,0._dp,routed_scalar_reset,accepted_start,good)
          if(.not.good) return
          scalar_reset=accepted_start
          correction=correction+closure_reset
          call average_down_reactive_eb_tree_3d(species,work,good,frame,held)
          if(.not.good) return
        endif
        call compute_reactive_eb_held_projection(0._dp,records(closing_record)%held_clock_end_change,0._dp, &
          scalar_reset,accepted_start,good)
        if(.not.good) return
        records(closing_record)%held_clock_end_change=accepted_start
        records(closing_record)%held_clock_events_closed=.true.
      enddo
    enddo
    ! The root exit belongs to the caller's still-pending parent closure.
    call average_down_reactive_eb_tree_field_3d(work,ir,good,minimum_depth=depth)
    if(.not.good) return
    if(present(heat_release)) then
      call average_down_reactive_eb_tree_field_3d(work,heat,good,minimum_depth=depth)
      if(.not.good) return
    endif
    if(.not.ieee_is_finite(correction).or.any(.not.ieee_is_finite(chemical))) return
    call allocate_held_mol_faces_3d(tree,nv,result);call allocate_held_mol_faces_3d(tree,1,scalar_result)
    do r=1,size(records)
      lev=records(r)%depth
      ! Parent closure uses adjacent levels, but the published physical
      ! family must retain every accepted descendant interval.
      do p=1,size(tree%nodes)
        if(tree%topology%nodes(p)%level/=lev) cycle
        result(p)%x=result(p)%x+family(r)%faces(p)%x
        result(p)%y=result(p)%y+family(r)%faces(p)%y
        result(p)%z=result(p)%z+family(r)%faces(p)%z
        scalar_result(p)%x=scalar_result(p)%x+family(r)%scalar_faces(p)%x
        scalar_result(p)%y=scalar_result(p)%y+family(r)%scalar_faces(p)%y
        scalar_result(p)%z=scalar_result(p)%z+family(r)%scalar_faces(p)%z
      enddo
    enddo
    tree=work;stored=held;reaction_source=ir;prefix=correction;chemical_change=chemical;history=records
    call move_alloc(result,flux);call move_alloc(scalar_result,scalar_flux)
    if(present(heat_release)) heat_release=heat
    ok=.true.
  end subroutine replay_reactive_eb_mol_held_group_3d

  ! Recursive SDC with interval-owned redistribution registers. Each register
  ! closes when the coarser source/recipient level finishes its interval.
  recursive subroutine advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,dt,iterations,chemistry, &
      transport_flags,force,rtol,atol,tree,extended_geometry,boundary_faces, &
      boundary_change,wall_change,chemical_change,source_change,minimum_theta,ok,reaction_source, &
      difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature,level_steps,redistribution, &
      target_volume_fraction,auxiliary_source,auxiliary_change,route_auxiliary, &
      evaluate_reaction,mol_mode,mol_lax_friedrichs,embedded_wall_temperature,embedded_noslip,frame, &
      state_redist_max_order,evaluate_spatial_diffusion,spatial_diffusion_halo_width, &
      constant_coefficients,heat_release,registered_flux,integrated_flux,internal_energy,internal_energy_flux, &
      energy_policy,energy_correction,held_interval_history,replay_interval_history,replay_spatial_weight, &
      replay_parent_registers,replay_parent_held_registers,provider_refused, &
      spatial_diffusion_is_source,auxiliary_has_internal_energy,internal_energy_pressure_work_change, &
      godunov_interval_history,godunov_replay_history,spatial_flux_weight, &
      godunov_parent_registers,godunov_parent_scalar_registers,godunov_routed_pending,godunov_replay_prefix, &
      godunov_recovery_history,godunov_recovery_ready,source_funding_context,native_replay_history, &
      native_spatial_weight,native_common_retries,native_replay_prefix,native_replay_preserve_weights, &
      separate_momentum,consistent_velocity,paired_pressure,connected_directions,observe_mol,observe_mol_sync, &
      observe_mol_reaction)
    procedure(eb_mol_sync_observer), optional :: observe_mol_sync
    logical :: observe_sync_enabled
    procedure(eb_mol_reaction_observer), optional :: observe_mol_reaction
    procedure(eb_mol_stage_observer), optional :: observe_mol
    logical, optional, intent(in) :: separate_momentum,consistent_velocity,paired_pressure,connected_directions
    logical :: selected_momentum,selected_velocity,selected_pressure,selected_connected
    ! Failure diagnostics; these do not commit a physical step or ledger.
    integer(int64), optional, intent(out) :: native_common_retries
    ! The full native receipt/open-face context is private with gas/H until
    ! outer success. Non-regular prepared funding still requires its native
    ! joint StateRedist stage path; this connects ordinary held MOL clocks.
    ! Rebuild a complete ordinary native MOL window from its clock anchor.
    ! External stages are frozen; gas/UEINT spatial operators, chemistry and
    ! receipts are rebuilt privately. Prefix continuation is explicitly opt-in.
    type(eb_mol_held_interval_3d), optional, intent(in) :: native_replay_history(:)
    real(dp), optional, intent(in) :: native_spatial_weight
    ! A prefix reuses only the supplied complete clocks; later clocks evaluate
    ! their sources natively. Publication still waits for the entire window.
    logical, optional, intent(in) :: native_replay_prefix
    logical, optional, intent(in) :: native_replay_preserve_weights
    logical :: replaying_native,native_prefix,preserve_native_weights
    integer :: native_cursor
    type(eb_held_funding_context_3d), optional, intent(inout) :: source_funding_context
    type(eb_held_funding_context_3d), allocatable :: funding_context_work
    logical, allocatable :: funding_local(:)
    integer :: funding_status
    type(eb_geometry_3d), allocatable :: funding_geometry(:)
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

    real(dp) :: godunov_weight,godunov_weight_low,godunov_weight_high
    type(eb_godunov_interval_3d), allocatable, optional, intent(inout) :: godunov_interval_history(:)
    type(eb_godunov_interval_3d), allocatable :: godunov_records(:)
    integer, allocatable :: godunov_active(:)
    integer :: godunov_cursor,godunov_capacity,godunov_level_calls,godunov_status
    logical :: godunov_capture,godunov_mol
    logical :: godunov_automatic_repair,godunov_provider_rejected,godunov_source_exhausted
    type :: godunov_root_trial
      type(eb_reactive_tree_3d) :: gas
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
    integer :: naux,auxiliary_alloc_status
    logical, optional, intent(in) :: spatial_diffusion_is_source
    logical :: external_diffusion_source
    type(eb_mol_held_interval_3d), optional, intent(in) :: replay_interval_history(:)
    real(dp), optional, intent(in) :: replay_spatial_weight
    type(state_redist_register_3d), allocatable, optional, intent(inout) :: replay_parent_registers(:)
    type(state_redist_register_3d), allocatable, optional, intent(inout) :: replay_parent_held_registers(:)
    logical, optional, intent(out) :: provider_refused
    type(eb_mol_held_interval_3d), allocatable, optional, intent(inout) :: held_interval_history(:)
    type(eb_mol_held_interval_3d), allocatable :: held_history(:)
    type :: held_replay_result
      type(eb_reactive_tree_3d) :: tree
      type(amr_tree_state_3d) :: stored,reaction_source
      type(eb_held_funding_context_3d), allocatable :: funding
      type(amr_tree_state_3d), allocatable :: heat
      type(eb_mol_held_interval_3d), allocatable :: records(:)
      type(eb_integrated_flux_3d), allocatable :: flux(:),scalar_flux(:)
      real(dp), allocatable :: chemical(:),wall(:),source(:),auxiliary(:)
      type(state_redist_register_3d), allocatable :: parent_registers(:),parent_held_registers(:)
      real(dp) :: prefix=0
    end type held_replay_result
    integer, allocatable :: deferred_record(:)
    logical :: held_replay_enabled,held_capture_enabled,held_replay_provider_refused,replaying_prepared
    integer :: replay_depth,replay_cursor
    real(dp) :: replay_time
    real(dp), allocatable :: replay_weight
    type(dual_energy_policy), optional, intent(in) :: energy_policy
    real(dp), optional, intent(inout) :: energy_correction
    type(dual_energy_policy) :: held_policy
    real(dp), allocatable :: correction_work
    real(dp) :: reset_change
    type(amr_tree_state_3d) :: reset_difference
    type(amr_tree_state_3d), optional, intent(inout) :: heat_release
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: internal_energy_flux(:)
    ! Time-integrated spatial face flux with pinned stage registration weights,
    ! before peer/coarse-fine synchronization. Separate from the physical flux
    ! used by the conservative state update and boundary inventory.
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: registered_flux(:)
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: integrated_flux(:)
    type(eb_integrated_flux_3d), allocatable :: registration(:)
    type(amr_tree_state_3d), allocatable :: heat
    real(dp), optional, intent(in) :: constant_coefficients(:)
    real(dp), optional, intent(in) :: embedded_wall_temperature
    logical, optional, intent(in) :: mol_mode,mol_lax_friedrichs,embedded_noslip
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    procedure(eb_sdc_spatial_diffusion_source), optional :: evaluate_spatial_diffusion
    integer, optional, intent(in) :: spatial_diffusion_halo_width
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: dt,force(3),rtol,atol
    integer, intent(in) :: iterations
    logical, intent(in) :: chemistry,transport_flags(4)
    type(eb_reactive_tree_3d), intent(inout) :: tree
    type(eb_geometry_3d), intent(in) :: extended_geometry(:)
    character(len=*), intent(in) :: boundary_faces(6)
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemical_change(:),source_change(:),minimum_theta
    integer, optional, intent(in) :: state_redist_max_order
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source
    real(dp), optional, intent(in) :: difmag,inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: wall_temperature(6)
    integer(int64), optional, intent(inout) :: level_steps(:)
    character(len=*), optional, intent(in) :: redistribution
    real(dp), optional, intent(in) :: target_volume_fraction
    procedure(eb_sdc_recursive_source), optional :: auxiliary_source
    real(dp), optional, intent(inout) :: auxiliary_change(:)
    real(dp), allocatable :: auxiliary(:)
    real(dp) :: source_time
    logical, optional, intent(in) :: route_auxiliary
    logical :: routing,old_auxiliary_only,has_auxiliary
    type(amr_tree_state_3d) :: off_old,off_new
    type(amr_tree_state_3d), allocatable :: auxiliary_pending(:)
    integer :: redist_order
    real(dp) :: redist_target
    character(len=16) :: method
    type :: interval_register_set
      type(state_redist_register_3d), allocatable :: fields(:)
    end type interval_register_set
    type(interval_register_set), allocatable :: interface_registers(:),held_interface_registers(:)
    logical :: prepared_redist,prepared_godunov
    integer :: godunov_padding
    type(eb_geometry_3d), allocatable :: register_geometry(:)
    type(amr_tree_state_3d), allocatable :: pending(:)
    type(eb_reactive_tree_3d) :: candidate,starts,ends
    type(amr_tree_state_3d) :: held_candidate,held_starts,held_ends
    type(amr_tree_state_3d), allocatable :: ir
    type(amr_tree_state_3d) :: source0,source1
    type(eb_integrated_flux_3d), allocatable :: total(:),synchronization_flux(:),held_total(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:),counts(:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: clock_start(:),clock_dt(:)
    real(dp) :: boundary(size(boundary_change)),wall(size(wall_change)),chemical(size(chemical_change))
    real(dp) :: external(size(source_change)),theta
    integer :: nv,node,n(3),levels,c
    logical :: valid
    if(present(native_common_retries)) native_common_retries=0_int64
    if(present(godunov_recovery_ready)) godunov_recovery_ready=.false.
    held_auxiliary=.false.
    if(present(auxiliary_has_internal_energy)) held_auxiliary=auxiliary_has_internal_energy
    external_diffusion_source=.false.
    if(present(spatial_diffusion_is_source)) external_diffusion_source=spatial_diffusion_is_source
    if(external_diffusion_source.and..not.present(evaluate_spatial_diffusion)) then
      ok=.false.;return
    endif
    observe_sync_enabled=present(observe_mol_sync)
    ok=.false.;nv=reactive_nvar(size(species))
    replaying_native=present(native_replay_history);native_cursor=0
    native_prefix=.false.
    if(present(native_replay_prefix)) native_prefix=native_replay_prefix
    preserve_native_weights=.false.
    if(present(native_replay_preserve_weights)) preserve_native_weights=native_replay_preserve_weights
    if(preserve_native_weights.and.(.not.replaying_native.or.present(native_spatial_weight))) return
    if(native_prefix.and..not.replaying_native) return
    if(present(native_spatial_weight)) then
      if(.not.replaying_native) return
      if(.not.ieee_is_finite(native_spatial_weight)) return
      if(native_spatial_weight<0.or.native_spatial_weight>1) return
    endif
    initial_clock_reset=0
    godunov_capture=present(godunov_interval_history).or.present(godunov_recovery_history)
    godunov_cursor=0;godunov_replay_cursor=0
    godunov_replay_depth=0;godunov_replay_base=0;godunov_parent=0
    godunov_replay_time=0;godunov_pressure_offset=0;godunov_minimum_offset=1
    godunov_provider_rejected=.false.;godunov_source_exhausted=.false.
    godunov_prefix=.false.
    if(present(godunov_replay_prefix)) godunov_prefix=godunov_replay_prefix
    if(godunov_prefix.and..not.present(godunov_replay_history)) return
    if(present(godunov_recovery_history).neqv.present(godunov_recovery_ready)) return
    godunov_weight=1
    if(present(spatial_flux_weight)) godunov_weight=spatial_flux_weight

    godunov_mol=.false.
    if(present(mol_mode)) godunov_mol=mol_mode
    selected_momentum=.false.;selected_velocity=.false.;selected_pressure=.false.;selected_connected=.false.
    if(present(separate_momentum)) selected_momentum=separate_momentum
    if(present(consistent_velocity)) selected_velocity=consistent_velocity
    if(present(paired_pressure)) selected_pressure=paired_pressure
    if(present(connected_directions)) selected_connected=connected_directions
    if(selected_velocity.and..not.selected_momentum) return
    if(selected_connected.and..not.selected_momentum) return
    if(selected_pressure.and..not.selected_velocity) return
    if(selected_pressure.and.present(frame)) return
    if(.not.ieee_is_finite(godunov_weight)) return
    if(godunov_weight<0.or.godunov_weight>1) return
    naux=eb_auxiliary_nvar(nv,auxiliary_has_internal_energy)
    allocate(auxiliary(naux),stat=auxiliary_alloc_status)
    if(auxiliary_alloc_status/=0) return
    if(held_auxiliary.and.(.not.present(internal_energy).or. &
      (.not.present(auxiliary_source).and..not.present(replay_interval_history)))) return
    if(present(provider_refused)) provider_refused=.false.
    has_auxiliary=present(auxiliary_source)
    replaying_prepared=present(replay_interval_history);replay_depth=0;replay_time=0;replay_cursor=0
    if(replaying_native.and.(replaying_prepared.or.godunov_capture.or.present(godunov_replay_history).or. &
      present(spatial_flux_weight))) return
    if(present(replay_spatial_weight).neqv.replaying_prepared) return
    if(replaying_prepared) then
      if(size(replay_interval_history)<1) return
      has_auxiliary=replay_interval_history(1)%has_auxiliary_source
      if(present(auxiliary_source).and..not.has_auxiliary) return
    endif
    if((present(replay_parent_registers).or.present(replay_parent_held_registers)).and..not.replaying_prepared) return
    held_policy=dual_energy_policy()
    if(present(energy_policy)) then
      if(.not.present(energy_correction)) return
      held_policy=energy_policy
    endif
    if(present(energy_correction)) then
      if(.not.present(internal_energy)) return
      if(.not.ieee_is_finite(energy_correction)) return
      allocate(correction_work,source=energy_correction)
    endif
    if(present(constant_coefficients)) then
      if(.not.valid_constant_transport_coefficients(constant_coefficients)) return
    endif
    redist_order=0
    if(present(state_redist_max_order)) redist_order=state_redist_max_order
    if(redist_order/=0.and.redist_order/=2.and.redist_order/=4) return
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    if(present(evaluate_reaction).and..not.chemistry) return
    routing=.false.
    if(present(route_auxiliary)) routing=route_auxiliary
    old_auxiliary_only=.false.
    if(present(mol_mode)) old_auxiliary_only=mol_mode.and.chemistry.and.iterations>1
    if(.not.valid_tree(tree,nv).or..not.ieee_is_finite(dt).or.dt<=0) return
    if(any([size(boundary_change),size(wall_change),size(chemical_change),size(source_change)]/=nv)) return
    if(routing.and..not.present(auxiliary_source)) return
    if(.not.replaying_prepared) then
      if(has_auxiliary.neqv.present(auxiliary_change)) return
    endif
    if(present(auxiliary_change)) then
      if(size(auxiliary_change)/=naux) return
    end if
    if(present(internal_energy_pressure_work_change).and..not.present(internal_energy)) return
    if(present(internal_energy).neqv.present(internal_energy_flux)) return
    if(present(internal_energy)) then
      if(has_auxiliary.and.((.not.chemistry.or.routing).and..not.held_auxiliary)) return
      if(chemistry.and..not.present(energy_correction)) return
      if(.not.internal_energy%is_valid(tree%topology)) return
      if(size(internal_energy%fields(1)%values,1)/=1) return
      do node=1,size(tree%nodes)
        if(.not.all(ieee_is_finite(internal_energy%fields(node)%values))) return
      enddo
    endif
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if(.not.valid) return
    levels=size(ratios)
    if(present(level_steps)) then
      if(size(level_steps)/=levels) return
    end if
    method='none'
    if(present(redistribution)) method=redistribution
    if(selected_momentum) then
      if(.not.godunov_mol.or..not.present(internal_energy).or.method/='none') return
    endif
    if(method/='none'.and.method/='flux_redist'.and.method/='state_redist') return
    prepared_redist=.false.
    if(present(mol_mode)) prepared_redist=mol_mode.and.method=='state_redist'
    if(present(internal_energy)) then
      if(method/='none'.and.method/='state_redist') return
    endif
    ! Grow9 or grow10 geometry selects the composed Godunov stage; grow6 callers retain
    ! their existing route until application geometry/callback migration.
    prepared_godunov=.false.
    if(method=='state_redist'.and..not.prepared_redist) then
      if(size(extended_geometry)==size(tree%nodes)) then
        prepared_godunov=.true.
        do node=1,size(tree%nodes)
          n=tree%topology%nodes(node)%shape
          if(.not.(all([extended_geometry(node)%nx,extended_geometry(node)%ny,extended_geometry(node)%nz]==n+18).or. &
            all([extended_geometry(node)%nx,extended_geometry(node)%ny,extended_geometry(node)%nz]==n+20))) &
            prepared_godunov=.false.
        enddo
      endif
    endif
    prepared_redist=prepared_redist.or.prepared_godunov
    if(present(internal_energy).and.method=='state_redist'.and..not.prepared_redist) return
    godunov_padding=9
    if(prepared_godunov.and.size(extended_geometry)>0) then
      godunov_padding=(extended_geometry(1)%nx-tree%topology%nodes(1)%shape(1))/2
    endif
    redist_target=.5_dp
    if(present(target_volume_fraction)) redist_target=target_volume_fraction
    if(.not.ieee_is_finite(redist_target)) return
    if(redist_target<=0.or.redist_target>1) return
    held_replay_enabled=.false.;held_replay_provider_refused=.false.
    if(present(mol_mode)) held_replay_enabled=mol_mode.and.chemistry.and.present(internal_energy).and. &
      (method=='none'.or.prepared_redist)
    if(replaying_prepared) then
      if(.not.held_replay_enabled.or..not.prepared_redist.or.prepared_godunov) return
      if(.not.ieee_is_finite(replay_spatial_weight).or.replay_spatial_weight<0.or.replay_spatial_weight>1) return
      call validate_prepared_clock_history(valid)
      if(.not.valid) return
      allocate(replay_weight,source=replay_spatial_weight)
      held_replay_enabled=.false.
    endif
    if(replaying_native) then
      if(.not.held_replay_enabled.or.method/='none'.or.routing) return
      call validate_native_window_history(valid)
      if(.not.valid) return
      held_replay_enabled=.false.
    endif
    held_capture_enabled=held_replay_enabled.or.replaying_prepared.or.replaying_native
    if(present(mol_mode).and.present(held_interval_history)) then
      held_capture_enabled=held_capture_enabled.or.(mol_mode.and.chemistry.and.present(internal_energy).and.prepared_redist)
    endif
    if(present(held_interval_history).and..not.held_capture_enabled) return
    if(held_capture_enabled) then
      allocate(held_history(0),deferred_record(0:levels-1));deferred_record=0
    endif
    if(godunov_capture) then
      if(godunov_mol.or..not.present(internal_energy).or.present(held_interval_history).or.replaying_prepared) return
      if(method/='none'.and..not.prepared_godunov) return
    endif
    if(present(godunov_replay_history).or.present(spatial_flux_weight)) then
      valid=.not.godunov_mol.and.present(internal_energy).and.(method=='none'.or.prepared_godunov)
      valid=valid.and..not.replaying_prepared.and..not.present(held_interval_history)
      if(.not.valid) return
    endif
    if(present(godunov_replay_history)) then
      call initialize_godunov_replay(valid)
      if(.not.valid) return
    else
      valid=.not.present(godunov_parent_registers).and..not.present(godunov_parent_scalar_registers).and. &
        .not.present(godunov_routed_pending)
      if(.not.valid) return
    endif
    if(present(godunov_replay_history)) then
      call validate_godunov_family(valid)
      if(.not.valid) return
    endif
    godunov_automatic_repair=.not.godunov_mol.and.chemistry.and.present(internal_energy).and. &
      .not.present(godunov_replay_history).and.(method=='none'.or.prepared_godunov)
    godunov_automatic_repair=godunov_automatic_repair.and..not.present(held_interval_history).and..not.replaying_prepared
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
      if(godunov_status/=0) return
      godunov_active=0
      if(godunov_replay_depth>0) godunov_active(godunov_replay_depth-1)=godunov_parent
    endif
    candidate=tree
    if(present(internal_energy)) held_candidate=internal_energy
    if(replaying_native) then
      candidate=native_replay_history(1)%clock_initial
      held_candidate=native_replay_history(1)%clock_internal_energy
      correction_work=native_replay_history(1)%clock_energy_correction
    endif
    if(replaying_prepared) then
      candidate=replay_interval_history(1)%clock_initial
      held_candidate=replay_interval_history(1)%clock_internal_energy
      correction_work=replay_interval_history(1)%clock_energy_correction
    else if(godunov_replay_depth==0.and..not.replaying_native) then
    if(allocated(correction_work)) then
      held_candidate=internal_energy
      call reset_held_eb_tree_family_3d(species,candidate,held_candidate,reset_difference,held_policy,-1, &
        reset_change,valid,frame, &
        internal_energy_change=initial_clock_reset,reset_level=merge(-2,-1,godunov_mol.and.chemistry.and.method=='none'))
      if(.not.valid) return
      correction_work=correction_work+reset_change
    endif
    if(.not.(godunov_mol.and.chemistry.and.method=='none'.and.allocated(correction_work))) then
    if(present(internal_energy)) then
      call average_down_reactive_eb_tree_3d(species,candidate,valid,frame,held_candidate)
    else
      call average_down_reactive_eb_tree_3d(species,candidate,valid,frame)
    endif
    if(.not.valid) return
    endif
    endif
    if(godunov_replay_depth>0) then
      associate(anchor=>godunov_replay_history(1))
      do node=1,size(tree%nodes)
        candidate%nodes(node)%state=anchor%clock_gas%fields(node)%values
        candidate%nodes(node)%temperature=anchor%clock_temperature%fields(node)%values(1,:,:,:)
      enddo
      call copy_godunov_field_3d(anchor%clock_internal_energy,held_candidate,valid)
      if(.not.valid) return
      if(allocated(correction_work)) correction_work=anchor%clock_energy_correction
      end associate
    endif
    if(present(source_funding_context)) then
      if(.not.present(mol_mode)) return
      if(.not.mol_mode.or..not.present(internal_energy).or..not.present(energy_correction)) return
      if(method/='none') return
      if(replaying_prepared.or.present(godunov_replay_history)) return
      allocate(funding_context_work,funding_local(size(tree%nodes)),funding_geometry(size(tree%nodes)), &
        stat=funding_status)
      if(funding_status/=0) return
      funding_local=.true.
      do node=1,size(tree%nodes)
        funding_geometry(node)=tree%nodes(node)%geometry
      enddo
      if(replaying_native) then
        call copy_eb_held_funding_context_3d(tree%topology,funding_geometry,nv,funding_local, &
          native_replay_history(1)%clock_source_funding_context,funding_context_work,valid)
      else if(allocated(source_funding_context%pending)) then
        call copy_eb_held_funding_context_3d(tree%topology,funding_geometry,nv,funding_local,source_funding_context, &
          funding_context_work,valid)
      else
        valid=.not.allocated(source_funding_context%open).and. &
          .not.allocated(source_funding_context%geometry).and..not.allocated(source_funding_context%topology%nodes)
        if(.not.valid) return
        call initialize_eb_held_funding_context_3d(tree%topology,funding_geometry,nv,funding_local, &
          funding_context_work,valid)
      endif
      if(.not.valid) return
      valid=valid_eb_held_funding_context_3d(tree%topology,funding_geometry,nv,funding_local, &
        funding_context_work,replay_depth)
      if(.not.valid) return
    endif
    starts=candidate;ends=candidate
    if(present(internal_energy)) then
      held_starts=held_candidate;held_ends=held_candidate
    endif
    allocate(clock_start(levels),clock_dt(levels),counts(levels),source0%fields(size(tree%nodes)))
    counts=0;clock_start=0;clock_dt=dt
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      allocate(source0%fields(node)%values(nv,n(1),n(2),n(3)),source=0._dp)
    end do
    source1=source0
    if(replaying_native) then
      associate(anchor=>native_replay_history(1))
      clock_start=anchor%clock_start;clock_dt=anchor%clock_interval
      starts=anchor%ancestor_start;ends=anchor%ancestor_end
      held_starts=anchor%ancestor_internal_energy_start;held_ends=anchor%ancestor_internal_energy_end
      source0=anchor%ancestor_source_start;source1=anchor%ancestor_source_end
      end associate
    endif
    if(godunov_replay_depth>0) then
      associate(anchor=>godunov_replay_history(1))
      clock_start=anchor%clock_start;clock_dt=anchor%clock_interval
      do node=1,size(tree%nodes)
        starts%nodes(node)%state=anchor%ancestor_start%fields(node)%values
        ends%nodes(node)%state=anchor%ancestor_end%fields(node)%values
        starts%nodes(node)%temperature=anchor%ancestor_temperature_start%fields(node)%values(1,:,:,:)
        ends%nodes(node)%temperature=anchor%ancestor_temperature_end%fields(node)%values(1,:,:,:)
      enddo
      call copy_godunov_field_3d(anchor%ancestor_internal_energy_start,held_starts,valid)
      if(.not.valid) return
      call copy_godunov_field_3d(anchor%ancestor_internal_energy_end,held_ends,valid)
      if(.not.valid) return
      call copy_godunov_field_3d(anchor%source_start,source0,valid)
      if(.not.valid) return
      call copy_godunov_field_3d(anchor%source_end,source1,valid)
      if(.not.valid) return
      end associate
    endif

    if(replaying_prepared.and.replay_depth>0) then
      starts=replay_interval_history(1)%ancestor_start;ends=replay_interval_history(1)%ancestor_end
      held_starts=replay_interval_history(1)%ancestor_internal_energy_start
      held_ends=replay_interval_history(1)%ancestor_internal_energy_end
      source0=replay_interval_history(1)%ancestor_source_start
      source1=replay_interval_history(1)%ancestor_source_end
    endif
    if(routing) then
      allocate(auxiliary_pending(0:levels-1))
      do c=0,levels-1
        auxiliary_pending(c)=source0
        if(held_auxiliary) then
          do node=1,size(tree%nodes)
            n=tree%topology%nodes(node)%shape
            deallocate(auxiliary_pending(c)%fields(node)%values)
            allocate(auxiliary_pending(c)%fields(node)%values(naux,n(1),n(2),n(3)),source=0._dp)
          enddo
        endif
      end do
    end if
    if(godunov_replay_depth>0.and.routing) then
      do c=0,levels-1
        call copy_godunov_field_3d(godunov_replay_history(1)%clock_routed_pending(c),auxiliary_pending(c),valid)
        if(.not.valid) return
      enddo
    endif
    if(method/='none'.and..not.prepared_redist) then
      allocate(pending(0:levels-1))
      do c=0,levels-1
        pending(c)=source0
      end do
    end if
    if(prepared_redist) then
      allocate(register_geometry(size(tree%nodes)))
      valid=.true.
      do node=1,size(tree%nodes)
        n=tree%topology%nodes(node)%shape
        if(prepared_godunov) then
          if(any([extended_geometry(node)%nx,extended_geometry(node)%ny,extended_geometry(node)%nz]/= &
            n+2*godunov_padding)) return
          call extract_eb_geometry_3d(extended_geometry(node),[1,1,1]*(godunov_padding-5), &
            n+godunov_padding+6,register_geometry(node),valid)
          if(.not.valid) exit
        else
          register_geometry(node)=extended_geometry(node)
        endif
      enddo
      if(.not.valid) return
      allocate(interface_registers(0:levels-2))
      if(present(internal_energy)) allocate(held_interface_registers(0:levels-2))
      do c=0,levels-2
        allocate(interface_registers(c)%fields(size(tree%nodes)))
        if(present(internal_energy)) allocate(held_interface_registers(c)%fields(size(tree%nodes)))
        do node=1,size(tree%nodes)
          call initialize_state_redist_tree_register_3d(tree%topology,node,nv,6, &
            interface_registers(c)%fields(node),valid,boundary_faces,scale=merge(1._dp,.5_dp,prepared_godunov))
          if(.not.valid) return
          if(present(internal_energy)) then
            call initialize_state_redist_tree_register_3d(tree%topology,node,1,6, &
              held_interface_registers(c)%fields(node),valid,boundary_faces,scale=merge(1._dp,.5_dp,prepared_godunov))
            if(.not.valid) return
          endif
        enddo
      enddo
    endif
    if(replaying_prepared.and.replay_depth>0) then
      interface_registers(replay_depth-1)%fields=replay_interval_history(1)%open_ancestor_registers
      held_interface_registers(replay_depth-1)%fields=replay_interval_history(1)%open_ancestor_held_registers
    endif
    if(godunov_replay_depth>0.and.prepared_redist) then
      associate(anchor=>godunov_replay_history(1))
      valid=.true.
      do node=1,size(tree%nodes)
        n=tree%topology%nodes(node)%shape
        godunov_expected_register=interface_registers(godunov_replay_depth-1)%fields(node)
        call accumulate_state_redist_register_3d(godunov_expected_register,anchor%open_ancestor_registers(node), &
          nv,n+12,valid)
        if(.not.valid) exit
        godunov_expected_register=held_interface_registers(godunov_replay_depth-1)%fields(node)
        call accumulate_state_redist_register_3d(godunov_expected_register,anchor%open_ancestor_scalar_registers(node), &
          1,n+12,valid)
        if(.not.valid) exit
      enddo
      if(.not.valid) return
      call copy_godunov_registers_3d(anchor%open_ancestor_registers,interface_registers(godunov_replay_depth-1)%fields,valid)
      if(.not.valid) return
      call copy_godunov_registers_3d(anchor%open_ancestor_scalar_registers, &
        held_interface_registers(godunov_replay_depth-1)%fields,valid)
      if(.not.valid) return
      end associate
    endif
    if(present(heat_release)) then
      if(.not.heat_release%is_valid(tree%topology)) return
      if(size(heat_release%fields(1)%values,1)/=1) return
      allocate(heat);heat=heat_release
    endif
    if(chemistry.or.present(reaction_source).or.godunov_capture.or.present(godunov_replay_history)) then
      allocate(ir)
      if(present(reaction_source)) then
        ir=reaction_source
      else
        allocate(ir%fields(size(tree%nodes)))
        do node=1,size(tree%nodes)
          n=tree%topology%nodes(node)%shape
          allocate(ir%fields(node)%values(nv-4,n(1),n(2),n(3)),source=0._dp)
        end do
      end if
    end if
    if(replaying_prepared) then
      ir=replay_interval_history(1)%clock_reaction_source
      if(present(heat_release)) heat=replay_interval_history(1)%clock_heat
    endif
    if(replaying_native) then
      ir=native_replay_history(1)%clock_reaction_source
      if(present(heat_release)) heat=native_replay_history(1)%clock_heat
    endif
    if(godunov_replay_depth>0) then
      call copy_godunov_field_3d(godunov_replay_history(1)%clock_reaction_source,ir,valid)
      if(.not.valid) return
      if(present(heat_release)) then
        call copy_godunov_field_3d(godunov_replay_history(1)%clock_heat,heat,valid)
        if(.not.valid) return
      endif
    endif
    auxiliary=0;source_time=0
    wall=0;chemical=0;external=0;theta=godunov_minimum_offset;pressure_work_change=0
    if(present(godunov_replay_history)) then
      replay_depth=godunov_replay_depth;replay_time=godunov_replay_time
    endif
    call advance_level(replay_depth,replay_time,dt,synchronization_flux,total,held_total,valid)
    if(present(provider_refused)) provider_refused=held_replay_provider_refused.or.godunov_provider_rejected
    if(.not.valid) then
      if(godunov_automatic_repair.and..not.godunov_provider_rejected.and..not.godunov_source_exhausted) &
        call repair_godunov_root()
      if(.not.ok.and.present(godunov_recovery_history)) call publish_godunov_recovery()
      return
    endif
    if(present(godunov_replay_history)) then
      if(godunov_replay_cursor/=godunov_replay_count) return
    endif
    if(replaying_prepared) then
      if(replay_cursor/=size(replay_interval_history)) return
    endif
    if(replaying_native) then
      if(native_cursor/=size(native_replay_history)) return
    endif
    if(allocated(correction_work)) then
      if(.not.(godunov_mol.and.chemistry.and.method=='none')) then
        call average_down_reactive_eb_tree_3d(species,candidate,valid,frame,held_candidate)
        if(.not.valid) return
      endif
      if(.not.ieee_is_finite(correction_work)) return
    endif
    if(prepared_redist) then
      ! All interval ledgers must be consumed before an outer step publishes.
      do c=replay_depth,levels-2
        do node=1,size(tree%nodes)
          if(allocated(interface_registers(c)%fields(node)%coarse_delta)) then
            if(any(interface_registers(c)%fields(node)%coarse_delta/=0)) return
          endif
          if(allocated(interface_registers(c)%fields(node)%fine_delta)) then
            if(any(interface_registers(c)%fields(node)%fine_delta/=0)) return
          endif
        enddo
      enddo
      if(present(internal_energy)) then
        do c=replay_depth,levels-2
          do node=1,size(tree%nodes)
            if(allocated(held_interface_registers(c)%fields(node)%coarse_delta)) then
              if(any(held_interface_registers(c)%fields(node)%coarse_delta/=0)) return
            endif
            if(allocated(held_interface_registers(c)%fields(node)%fine_delta)) then
              if(any(held_interface_registers(c)%fields(node)%fine_delta/=0)) return
            endif
          enddo
        enddo
      endif
    endif
    n=tree%topology%nodes(1)%shape
    associate(g=>tree%nodes(1)%geometry,f=>total(1))
    do c=1,nv
      boundary(c)=g%dy*g%dz*sum(g%x_face_fraction(0,:,:)*f%x(c,0,:,:)- &
        g%x_face_fraction(n(1),:,:)*f%x(c,n(1),:,:))+ &
        g%dx*g%dz*sum(g%y_face_fraction(:,0,:)*f%y(c,:,0,:)- &
        g%y_face_fraction(:,n(2),:)*f%y(c,:,n(2),:))+ &
        g%dx*g%dy*sum(g%z_face_fraction(:,:,0)*f%z(c,:,:,0)- &
        g%z_face_fraction(:,:,n(3))*f%z(c,:,:,n(3)))
    end do
    end associate
    if(.not.all(ieee_is_finite(boundary)).or..not.all(ieee_is_finite(wall)).or. &
      .not.all(ieee_is_finite(chemical)).or..not.all(ieee_is_finite(external))) return
    if(godunov_capture) then
      call rebuild_godunov_held_clock_prefixes_3d(godunov_records,valid)
      if(.not.valid) return
    endif
    if(held_capture_enabled) then
      call rebuild_serial_mol_held_clock_prefixes_3d(held_history,valid)
      if(.not.valid) return
    endif
    call prepare_godunov_parent_outputs(valid)
    if(.not.valid) return
    if(replaying_prepared) then
      if(replay_depth>0) then
        replay_parent_registers=interface_registers(replay_depth-1)%fields
        replay_parent_held_registers=held_interface_registers(replay_depth-1)%fields
      else
        if(present(replay_parent_registers)) then
          if(allocated(replay_parent_registers)) deallocate(replay_parent_registers)
        endif
        if(present(replay_parent_held_registers)) then
          if(allocated(replay_parent_held_registers)) deallocate(replay_parent_held_registers)
        endif
      endif
    endif
    if(present(held_interval_history)) held_interval_history=held_history
    if(godunov_capture) then
      if(godunov_cursor/=size(godunov_records)) return
      if(present(godunov_interval_history)) call move_alloc(godunov_records,godunov_interval_history)
    endif
    if(present(godunov_parent_registers)) call move_alloc(godunov_parent_work,godunov_parent_registers)
    if(present(godunov_parent_scalar_registers)) &
      call move_alloc(godunov_parent_scalar_work,godunov_parent_scalar_registers)
    if(present(godunov_routed_pending)) call move_alloc(godunov_pending_work,godunov_routed_pending)
    if(present(energy_correction)) energy_correction=correction_work
    if(present(internal_energy_pressure_work_change)) internal_energy_pressure_work_change=pressure_work_change
    tree=candidate;boundary_change=boundary;wall_change=wall;chemical_change=chemical;source_change=external
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
  contains
    subroutine apply_native_funding_receipt(receipt,good)
      type(amr_tree_state_3d), intent(in) :: receipt
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: gas_work,held_work
      integer :: p,allocation_status
      good=.false.
      allocate(gas_work%fields(size(tree%nodes)),stat=allocation_status)
      if(allocation_status/=0) return
      do p=1,size(tree%nodes)
        allocate(gas_work%fields(p)%values,source=candidate%nodes(p)%state,stat=allocation_status)
        if(allocation_status/=0) return
      enddo
      call copy_godunov_field_3d(held_candidate,held_work,good);if(.not.good) return
      call apply_eb_held_auxiliary_increment_3d(gas_work,held_work,receipt,good);if(.not.good) return
      do p=1,size(tree%nodes)
        call move_alloc(gas_work%fields(p)%values,candidate%nodes(p)%state)
      enddo
      call move_alloc(held_work%fields,held_candidate%fields)
    end subroutine apply_native_funding_receipt

    subroutine initialize_godunov_replay(good)
      logical, intent(out) :: good
      logical :: local_ok
      integer :: lev,calls
      local_ok=size(godunov_replay_history)>0.and.lbound(godunov_replay_history,1)==1
      good=local_ok;if(.not.good) return
      associate(anchor=>godunov_replay_history(1))
      local_ok=anchor%depth>=0.and.anchor%depth<levels.and. &
        anchor%subtree_last>=size(godunov_replay_history).and.ieee_is_finite(anchor%time).and. &
        ieee_is_finite(anchor%interval)
      good=local_ok;if(.not.good) return
      godunov_replay_depth=anchor%depth;godunov_replay_time=anchor%time
      godunov_replay_count=1;calls=1
      do lev=godunov_replay_depth+2,levels
        local_ok=calls<=huge(calls)/ratios(lev)
      good=local_ok;if(.not.good) return
        calls=calls*ratios(lev)
        local_ok=godunov_replay_count<=huge(calls)-calls
      good=local_ok;if(.not.good) return
        godunov_replay_count=godunov_replay_count+calls
      enddo
      local_ok=size(godunov_replay_history)<=godunov_replay_count
      if(.not.godunov_prefix) local_ok=size(godunov_replay_history)==godunov_replay_count
      good=local_ok;if(.not.good) return
      godunov_replay_base=anchor%subtree_last-godunov_replay_count;godunov_parent=anchor%parent_record
      local_ok=anchor%interval==dt.and.godunov_replay_time>=0
      if(godunov_replay_depth==0) then
        local_ok=local_ok.and.godunov_replay_base==0.and.godunov_parent==0.and.godunov_replay_time==0
      else
        local_ok=local_ok.and.godunov_replay_base>=1.and.godunov_parent>=1.and.godunov_parent<=godunov_replay_base
      endif
      good=local_ok;if(.not.good) return
      if(godunov_replay_depth>0) then
        local_ok=valid_godunov_clock_anchor_3d(anchor,tree%topology,nv,naux,levels)
      good=local_ok;if(.not.good) return
        local_ok=.true.
        if(prepared_redist) local_ok=present(godunov_parent_registers).and.present(godunov_parent_scalar_registers)
        if(routing) local_ok=local_ok.and.present(godunov_routed_pending)
        if(present(heat_release)) local_ok=local_ok.and.allocated(anchor%clock_heat%fields)
      good=local_ok;if(.not.good) return
        local_ok=ieee_is_finite(held_policy%eta2)
      good=local_ok;if(.not.good) return
        local_ok=anchor%policy%eta2==held_policy%eta2.and. &
          (anchor%policy%allow_small_energy.eqv.held_policy%allow_small_energy).and. &
          (anchor%policy%allow_negative_energy.eqv.held_policy%allow_negative_energy).and. &
          (anchor%policy%update_total_from_internal.eqv.held_policy%update_total_from_internal)
      good=local_ok;if(.not.good) return
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
      good=local_ok;if(.not.good) return
        call copy_godunov_registers_3d(held_interface_registers(godunov_replay_depth-1)%fields, &
          godunov_parent_scalar_work,local_ok)
      good=local_ok;if(.not.good) return
      endif
      if(present(godunov_routed_pending).and.routing) then
        allocate(godunov_pending_work(0:levels-1),stat=allocation_status)
        local_ok=allocation_status==0
      good=local_ok;if(.not.good) return
        do bucket=0,levels-1
          local_ok=.true.
          do p=1,size(tree%nodes)
            local_ok=local_ok.and.all(ieee_is_finite(auxiliary_pending(bucket)%fields(p)%values))
          enddo
      good=local_ok;if(.not.good) return
          call copy_godunov_field_3d(auxiliary_pending(bucket),godunov_pending_work(bucket),local_ok)
      good=local_ok;if(.not.good) return
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
      good=local_ok;if(.not.good) return
      cursor=first-1
      call seal_godunov_prefix_clock(godunov_records(first)%depth,godunov_records(first)%parent_record, &
        godunov_records(first)%time,godunov_records(first)%interval,cursor,last,good)
      if(.not.good) return
      local_ok=last<=cursor
      good=local_ok;if(.not.good) return
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
      good=local_ok;if(.not.good) return
      cursor=cursor+1;here=cursor
      if(here<=last) then
        associate(saved=>godunov_records(here))
        local_ok=valid_godunov_source_cache_3d(saved,tree%topology,nv,naux,allow_partial=.true.)
      good=local_ok;if(.not.good) return
        local_ok=saved%accepted_source_stages==iterations+1.or.here==last
      good=local_ok;if(.not.good) return
        tolerance=32*epsilon(1._dp)*max(abs(time+duration),tiny(1._dp))
        local_ok=saved%depth==depth.and.saved%parent_record==parent.and. &
          ieee_is_finite(saved%time).and.ieee_is_finite(saved%interval)
      good=local_ok;if(.not.good) return
        local_ok=abs(saved%time-time)<=tolerance.and.abs(saved%interval-duration)<=tolerance
      good=local_ok;if(.not.good) return
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
      if(godunov_provider_rejected.or.godunov_source_exhausted.or.held_replay_provider_refused) return
      local_ok=allocated(godunov_records).and.godunov_cursor>0
      ready=local_ok;if(.not.ready) return
      local_ok=godunov_cursor<=size(godunov_records)
      ready=local_ok;if(.not.ready) return
      call seal_godunov_group_prefix(1,godunov_cursor,ready)
      if(.not.ready) return
      allocate(work(godunov_cursor),stat=allocation_status)
      local_ok=allocation_status==0
      ready=local_ok;if(.not.ready) return
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
      good=local_ok;if(.not.good) return
      known=last-first+1
      recovered_last=first-1+size(result%recovery)
      local_ok=recovered_last<=size(godunov_records)
      good=local_ok;if(.not.good) return
      ! Keep original accepted packets and their original incoming clock.
      ! Add only genuinely new accepted tables from the quarantined trial.
      do r=1,size(result%recovery)
        local_ok=valid_godunov_source_cache_3d(result%recovery(r),tree%topology,nv,naux,allow_partial=.true.)
      good=local_ok;if(.not.good) return
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
      if(godunov_provider_rejected.or.godunov_source_exhausted.or.held_replay_provider_refused.or. &
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
      good=local_ok;if(.not.good) return
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
      good=local_ok;if(.not.good) return
      ! Prepare physical and outer registration families before adoption.
      call copy_godunov_fluxes_3d(best%faces,replacement_physical,local_ok)
      good=local_ok;if(.not.good) return
      if(present(registered_flux)) then
        call copy_godunov_fluxes_3d(registration,replacement_registration,local_ok)
      good=local_ok;if(.not.good) return
        do p=1,size(tree%nodes)
          replacement_registration(p)%x=0;replacement_registration(p)%y=0;replacement_registration(p)%z=0
        enddo
        do r=1,replacement_last
          if(r>=first) then
            lev=best%records(r-first+1)%depth
          else
            lev=godunov_records(r)%depth
          endif
          do p=1,size(tree%nodes)
            if(tree%topology%nodes(p)%level/=lev) cycle
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
          local_ok=local_ok.and.all(ieee_is_finite(replacement_registration(p)%x)).and. &
            all(ieee_is_finite(replacement_registration(p)%y)).and.all(ieee_is_finite(replacement_registration(p)%z))
        enddo
      good=local_ok;if(.not.good) return
      endif
      candidate=best%gas
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
        lev=tree%topology%nodes(p)%level
        if(lev<depth) cycle
        latest=0
        do r=1,size(best%records)
          if(best%records(r)%depth==lev) latest=r
        enddo
        if(latest==0) cycle
        associate(saved=>best%records(latest))
        starts%nodes(p)%state=saved%initial_gas%fields(p)%values
        ends%nodes(p)%state=saved%final_gas%fields(p)%values
        starts%nodes(p)%temperature=saved%initial_temperature%fields(p)%values(1,:,:,:)
        ends%nodes(p)%temperature=saved%final_temperature%fields(p)%values(1,:,:,:)
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
      valid=local_ok;if(.not.valid) return
      do r=1,size(godunov_records)
        local_ok=valid_godunov_source_cache_3d(godunov_records(r),tree%topology,nv,naux)
      valid=local_ok;if(.not.valid) return
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
      tree=best%gas
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
      good=local_ok;if(.not.good) return
      prefix_replay=size(frozen)<frozen(1)%subtree_last-first+1.or.any(frozen%source_prefix_only)
      result%gas=tree
      result%stored=internal_energy;result%prefix=energy_correction
      result%reaction=godunov_records(first)%clock_reaction_source
      if(present(reaction_source)) result%reaction=reaction_source
      allocate(result%boundary(nv),result%wall(nv),result%chemical(nv),result%external(nv), &
        result%counts(levels),stat=allocation_status)
      local_ok=allocation_status==0
      if(.not.local_ok) return
      if(has_auxiliary) then
        allocate(result%auxiliary(naux),stat=allocation_status)
        local_ok=allocation_status==0
        if(.not.local_ok) return
      endif
      if(present(heat_release)) then
        allocate(result%heat,source=heat_release,stat=allocation_status)
        local_ok=allocation_status==0
        if(.not.local_ok) return
      endif
      call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport, &
        godunov_records(first)%interval,iterations,chemistry, &
        transport_flags,force,rtol,atol,result%gas,extended_geometry,boundary_faces, &
        result%boundary,result%wall,result%chemical,result%external,result%limiter,good, &
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
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        observe_mol=observe_mol,observe_mol_sync=observe_mol_sync,observe_mol_reaction=observe_mol_reaction)
    end subroutine build_godunov_group_trial
    subroutine validate_godunov_family(good)
      logical, intent(out) :: good
      integer :: cursor
      logical :: local_ok
      cursor=0
      local_ok=lbound(godunov_replay_history,1)==1.and.size(godunov_replay_history)>0
      good=local_ok;if(.not.good) return
      call validate_godunov_clock(godunov_replay_depth,godunov_parent,godunov_replay_time,dt,cursor,good)
      if(.not.good) return
      local_ok=cursor==godunov_replay_count
      good=local_ok;if(.not.good) return
    end subroutine validate_godunov_family
    recursive subroutine validate_godunov_clock(depth,parent,time,duration,cursor,good)
      integer, intent(in) :: depth,parent
      real(dp), intent(in) :: time,duration
      integer, intent(inout) :: cursor
      logical, intent(out) :: good
      integer :: here,step
      real(dp) :: elapsed,child_dt,tolerance
      type(eb_geometry_3d), allocatable :: funding_geometry(:)
      integer :: funding_node,funding_status
      logical :: local_ok
      local_ok=cursor<godunov_replay_count
      good=local_ok;if(.not.good) return
      cursor=cursor+1;here=cursor
      if(here<=size(godunov_replay_history)) then
      associate(saved=>godunov_replay_history(here))
      local_ok=ieee_is_finite(saved%time).and.ieee_is_finite(saved%interval)
      good=local_ok;if(.not.good) return
      tolerance=32*epsilon(1._dp)*max(abs(time+duration),tiny(1._dp))
      local_ok=saved%depth==depth.and.saved%parent_record==parent.and. &
        abs(saved%time-time)<=tolerance.and.abs(saved%interval-duration)<=tolerance.and. &
        saved%iterations==iterations.and.(saved%chemistry.eqv.chemistry).and. &
        all(saved%transport_flags.eqv.transport_flags).and.all(ieee_is_finite(saved%force))
      good=local_ok;if(.not.good) return
      local_ok=all(saved%force==force).and.(saved%prepared_redistribution.eqv.prepared_redist).and. &
        (saved%has_auxiliary.eqv.has_auxiliary).and.(saved%has_spatial_provider.eqv.present(evaluate_spatial_diffusion)).and. &
        (saved%spatial_diffusion_is_source.eqv.external_diffusion_source).and. &
        (saved%auxiliary_has_internal_energy.eqv.held_auxiliary).and.(saved%routed_auxiliary.eqv.routing)
      good=local_ok;if(.not.good) return
      local_ok=.not.saved%has_valid_diffusion_provider
      good=local_ok;if(.not.good) return
      local_ok=saved%accepted_source_stages==iterations+1.or.here==size(godunov_replay_history)
      good=local_ok;if(.not.good) return
      if(allocated(saved%clock_source_funding_context%pending)) then
        allocate(funding_geometry(size(tree%nodes)),stat=funding_status)
        if(funding_status/=0) return
        do funding_node=1,size(tree%nodes)
          funding_geometry(funding_node)=tree%nodes(funding_node)%geometry
        enddo
        local_ok=valid_eb_held_funding_context_3d(tree%topology,funding_geometry,nv, &
          [( .true.,funding_node=1,size(tree%nodes))],saved%clock_source_funding_context,depth)
        good=local_ok;if(.not.good) return
      endif
      local_ok=valid_godunov_source_cache_3d(saved,tree%topology,nv,naux,allow_partial=godunov_prefix)
      good=local_ok;if(.not.good) return
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
      good=local_ok;if(.not.good) return
      endif
      good=.true.
    end subroutine validate_godunov_clock
    subroutine capture_godunov_clock(level,saved,good)
      integer, intent(in) :: level
      type(eb_godunov_interval_3d), intent(inout) :: saved
      logical, intent(out) :: good
      logical :: local_ok
      integer :: bucket
      call pack_godunov_tree_3d(candidate,saved%clock_gas,saved%clock_temperature,local_ok)
      good=local_ok;if(.not.good) return
      call copy_godunov_field_3d(ir,saved%clock_reaction_source,local_ok)
      good=local_ok;if(.not.good) return
      call copy_godunov_field_3d(held_candidate,saved%clock_internal_energy,local_ok)
      good=local_ok;if(.not.good) return
      if(present(heat_release)) then
      call copy_godunov_field_3d(heat,saved%clock_heat,local_ok)
      good=local_ok;if(.not.good) return
      endif
      saved%clock_energy_correction=0
      if(allocated(correction_work)) saved%clock_energy_correction=correction_work
      saved%clock_pressure_work_change=godunov_pressure_offset+pressure_work_change
      saved%clock_minimum_theta=theta
      saved%routed_auxiliary=routing
      if(routing) then
        allocate(saved%clock_routed_pending(0:levels-1),stat=godunov_status)
        local_ok=godunov_status==0
      good=local_ok;if(.not.good) return
        do bucket=0,levels-1
          call copy_godunov_field_3d(auxiliary_pending(bucket),saved%clock_routed_pending(bucket),local_ok)
      good=local_ok;if(.not.good) return
        enddo
      endif
      allocate(saved%clock_start(levels),saved%clock_interval(levels),stat=godunov_status)
      local_ok=godunov_status==0
      good=local_ok;if(.not.good) return
      saved%clock_start=clock_start;saved%clock_interval=clock_dt
      if(prepared_redist.and.level>0) then
      call copy_godunov_registers_3d(interface_registers(level-1)%fields,saved%open_ancestor_registers,local_ok)
      good=local_ok;if(.not.good) return
      call copy_godunov_registers_3d(held_interface_registers(level-1)%fields,saved%open_ancestor_scalar_registers,local_ok)
      good=local_ok;if(.not.good) return
      endif
      good=.true.
    end subroutine capture_godunov_clock
    subroutine close_godunov_clock(saved,faces,scalar_faces,good)
      type(eb_godunov_interval_3d), intent(inout) :: saved
      type(eb_integrated_flux_3d), allocatable, intent(in) :: faces(:),scalar_faces(:)
      logical, intent(out) :: good
      logical :: local_ok
      call pack_godunov_tree_3d(candidate,saved%closed_gas,saved%closed_temperature,local_ok)
      good=local_ok;if(.not.good) return
      call copy_godunov_field_3d(held_candidate,saved%closed_internal_energy,local_ok)
      good=local_ok;if(.not.good) return
      call copy_godunov_field_3d(ir,saved%closed_reaction_source,local_ok)
      good=local_ok;if(.not.good) return
      call copy_godunov_fluxes_3d(faces,saved%closed_faces,local_ok)
      good=local_ok;if(.not.good) return
      call copy_godunov_fluxes_3d(scalar_faces,saved%closed_scalar_faces,local_ok)
      good=local_ok;if(.not.good) return
      saved%subtree_last=godunov_cursor+godunov_replay_base;saved%closed_energy_correction=0
      if(allocated(correction_work)) saved%closed_energy_correction=correction_work
      saved%closed_pressure_work_change=godunov_pressure_offset+pressure_work_change
      good=.true.
    end subroutine close_godunov_clock
    subroutine validate_native_window_history(good)
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: gas_rate,held_rate
      type(eb_mol_held_interval_3d) :: snapshot
      real(dp) :: cursor(0:levels-1),ending,tolerance,expected_dt,expected_start,expected_end
      real(dp) :: parent_start(0:levels-1),parent_dt(0:levels-1)
      integer :: seen(0:levels-1),r,lev,c,p
      logical :: reconstructed
      good=.false.
      if(size(native_replay_history)<1) return
      associate(anchor=>native_replay_history(1))
      replay_depth=anchor%depth;replay_time=anchor%time
      if(.not.valid_native_mol_clock_profile_3d(anchor%clock_start,anchor%clock_interval, &
        replay_depth,replay_time,dt,levels)) return
      if(anchor%interval/=dt.or.(replay_depth==0.and.replay_time/=0)) return
      if(replay_depth>0.and..not.present(source_funding_context)) return
      parent_start=anchor%clock_start;parent_dt=anchor%clock_interval
      end associate
      ending=replay_time+dt
      if(.not.ieee_is_finite(ending).or.ending<=replay_time) return
      cursor=replay_time;seen=0
      tolerance=32*epsilon(1._dp)*max(abs(ending),tiny(1._dp))
      expected_dt=dt
      do c=replay_depth-1,0,-1
        if(expected_dt>huge(1._dp)/real(ratios(c+2),dp)) return
        expected_dt=expected_dt*real(ratios(c+2),dp)
        if(abs(parent_dt(c)-expected_dt)>tolerance) return
      enddo
      do r=1,size(native_replay_history)
        associate(saved=>native_replay_history(r))
        if(saved%separate_momentum.neqv.selected_momentum) return
        if(saved%consistent_velocity.neqv.selected_velocity) return
        if(saved%paired_pressure.neqv.selected_pressure) return
        if(saved%connected_directions.neqv.selected_connected) return
        if(.not.ieee_is_finite(saved%native_spatial_weight)) return
        if(saved%native_spatial_weight<0.or.saved%native_spatial_weight>1) return
        if(.not.valid_held_mol_record_3d(tree,saved,nv,frame)) return
        lev=saved%depth
        if(lev<replay_depth) return
        if(.not.valid_native_mol_clock_profile_3d(saved%clock_start,saved%clock_interval, &
          lev,saved%time,saved%interval,levels)) return
        if(saved%native_iterations/=iterations.or.any(saved%native_force/=force)) return
        if(saved%prepared_redistribution.or.saved%routed_auxiliary) return
        if(saved%has_auxiliary_source.neqv.has_auxiliary) return
        if(saved%has_source_funding.neqv.present(source_funding_context)) return
        if(saved%spatial_diffusion_is_source.neqv.external_diffusion_source) return
        if(saved%auxiliary_has_internal_energy.neqv.held_auxiliary) return
        if(any(saved%native_providers.neqv.[present(evaluate_spatial_diffusion),.false.])) return
        if(any(saved%native_providers).and..not.external_diffusion_source) return
        if(.not.valid_mol_native_provider_cache_3d(tree,saved,nv)) return
        if(saved%policy%eta2/=held_policy%eta2.or. &
          (saved%policy%allow_small_energy.neqv.held_policy%allow_small_energy).or. &
          (saved%policy%allow_negative_energy.neqv.held_policy%allow_negative_energy).or. &
          (saved%policy%update_total_from_internal.neqv.held_policy%update_total_from_internal)) return
        if(.not.allocated(saved%fractions_start).or..not.allocated(saved%fractions_end)) return
        if(size(saved%fractions_start)/=lev.or.size(saved%fractions_end)/=lev) return
        if(.not.all(ieee_is_finite(saved%fractions_start)).or. &
          .not.all(ieee_is_finite(saved%fractions_end))) return
        if(any(saved%fractions_start<0).or.any(saved%fractions_end>1).or. &
          any(saved%fractions_end<=saved%fractions_start)) return
        if(.not.valid_tree(saved%ancestor_start,nv).or..not.valid_tree(saved%ancestor_end,nv)) return
        snapshot=saved;snapshot%initial=saved%ancestor_start
        if(.not.valid_held_mol_record_3d(tree,snapshot,nv,frame)) return
        snapshot%initial=saved%ancestor_end
        if(.not.valid_held_mol_record_3d(tree,snapshot,nv,frame)) return
        if(.not.valid_held_mol_field_3d(tree,saved%ancestor_internal_energy_start,1).or. &
          .not.valid_held_mol_field_3d(tree,saved%ancestor_internal_energy_end,1).or. &
          .not.valid_held_mol_field_3d(tree,saved%ancestor_source_start,nv).or. &
          .not.valid_held_mol_field_3d(tree,saved%ancestor_source_end,nv)) return
        if(abs(saved%time-cursor(lev))>tolerance) return
        if(r>1) then
          if(lev>native_replay_history(r-1)%depth+1) return
        endif
        expected_dt=dt
        do c=replay_depth+2,lev+1
          expected_dt=expected_dt/real(ratios(c),dp)
        enddo
        if(abs(saved%interval-expected_dt)>tolerance) return
        do c=0,lev-1
          if(saved%clock_start(c+1)/=parent_start(c).or.saved%clock_interval(c+1)/=parent_dt(c)) return
          expected_start=(saved%time-parent_start(c))/parent_dt(c)
          expected_end=(saved%time+saved%interval-parent_start(c))/parent_dt(c)
          if(abs(saved%fractions_start(c+1)-expected_start)>32*epsilon(1._dp).or. &
            abs(saved%fractions_end(c+1)-expected_end)>32*epsilon(1._dp)) return
        enddo
        do p=1,size(tree%nodes)
          if(tree%topology%nodes(p)%level>=replay_depth) cycle
          associate(anchor=>native_replay_history(1))
          if(any(saved%ancestor_start%nodes(p)%state/=anchor%ancestor_start%nodes(p)%state).or. &
            any(saved%ancestor_end%nodes(p)%state/=anchor%ancestor_end%nodes(p)%state).or. &
            any(saved%ancestor_start%nodes(p)%temperature/=anchor%ancestor_start%nodes(p)%temperature).or. &
            any(saved%ancestor_end%nodes(p)%temperature/=anchor%ancestor_end%nodes(p)%temperature)) return
          if(any(saved%ancestor_internal_energy_start%fields(p)%values/= &
            anchor%ancestor_internal_energy_start%fields(p)%values).or. &
            any(saved%ancestor_internal_energy_end%fields(p)%values/= &
            anchor%ancestor_internal_energy_end%fields(p)%values).or. &
            any(saved%ancestor_source_start%fields(p)%values/=anchor%ancestor_source_start%fields(p)%values).or. &
            any(saved%ancestor_source_end%fields(p)%values/=anchor%ancestor_source_end%fields(p)%values)) return
          end associate
        enddo
        parent_start(lev)=saved%time;parent_dt(lev)=saved%interval
        cursor(lev)=saved%time+saved%interval;seen(lev)=seen(lev)+1
        if(cursor(lev)>ending+tolerance) return
        call reconstruct_reactive_eb_mol_held_interval_3d(species,tree,saved,saved%spatial_flux, &
          saved%internal_energy_flux,gas_rate,held_rate,reconstructed,frame)
        if(.not.reconstructed) return
        end associate
      enddo
      if(.not.native_prefix.and.(any(seen(replay_depth:)<1).or.any(abs(cursor(replay_depth:)-ending)>tolerance))) then
        good=.false.;return
      endif
      associate(anchor=>native_replay_history(1))
      snapshot=anchor;snapshot%initial=anchor%clock_initial
      good=valid_held_mol_record_3d(tree,snapshot,nv,frame).and. &
        valid_held_mol_field_3d(tree,anchor%clock_internal_energy,1).and. &
        valid_held_mol_field_3d(tree,anchor%clock_reaction_source,nv-4)
      if(present(heat_release)) good=good.and.valid_held_mol_field_3d(tree,anchor%clock_heat,1)
      end associate
    end subroutine validate_native_window_history

    subroutine validate_prepared_clock_history(good)
      logical, intent(out) :: good
      type(eb_mol_held_interval_3d) :: identity
      type(state_redist_register_3d) :: expected_register
      real(dp) :: cursor(0:levels-1),ending,tolerance,parent_dt,parent_start,expected_start,expected_end
      integer :: r,lev,p,c,d(3)
      integer :: seen(0:levels-1)
      logical :: valid
      good=.false.
      if(size(replay_interval_history)<1) return
      replay_depth=replay_interval_history(1)%depth;replay_time=replay_interval_history(1)%time
      if(replay_depth<0.or.replay_depth>=levels.or.replay_interval_history(1)%interval/=dt) return
      if(replay_depth>0) then
        if(.not.present(replay_parent_registers).or..not.present(replay_parent_held_registers)) return
      endif
      ending=replay_time+dt
      if(.not.ieee_is_finite(ending).or.ending<=replay_time) return
      tolerance=32*epsilon(1._dp)*max(abs(ending),tiny(1._dp));cursor=replay_time;seen=0
      do r=1,size(replay_interval_history)
        associate(saved=>replay_interval_history(r))
        lev=saved%depth
        if(lev<replay_depth.or.lev>=levels) return
        if(saved%spatial_diffusion_is_source.neqv.external_diffusion_source) return
        if(saved%auxiliary_has_internal_energy.neqv.held_auxiliary) return
        if(saved%routed_auxiliary.neqv.routing) return
        if(.not.valid_held_prepared_mol_record_3d(tree,saved,nv,lev,saved%interval,iterations, &
          redist_target,redist_order,held_policy,frame,has_auxiliary)) return
        if(abs(saved%time-cursor(lev))>tolerance) return
        cursor(lev)=saved%time+saved%interval;seen(lev)=seen(lev)+1
        if(.not.ieee_is_finite(cursor(lev)).or.cursor(lev)>ending+tolerance) return
        if(r>1) then
          if(lev>replay_interval_history(r-1)%depth+1) return
        endif
        if(.not.allocated(saved%fractions_start).or..not.allocated(saved%fractions_end)) return
        if(size(saved%fractions_start)/=lev.or.size(saved%fractions_end)/=lev) return
        if(any(.not.ieee_is_finite(saved%fractions_start)).or.any(.not.ieee_is_finite(saved%fractions_end))) return
        if(any(saved%fractions_start<0).or.any(saved%fractions_end>1).or. &
          any(saved%fractions_end<=saved%fractions_start)) return
        if(.not.allocated(saved%wall_change).or..not.allocated(saved%source_change).or. &
          .not.allocated(saved%auxiliary_change)) return
        if(size(saved%wall_change)/=nv.or.size(saved%source_change)/=nv.or.size(saved%auxiliary_change)/=naux) return
        if(any(.not.ieee_is_finite(saved%wall_change)).or.any(.not.ieee_is_finite(saved%source_change)).or. &
          any(.not.ieee_is_finite(saved%auxiliary_change))) return
        end associate
      enddo
      if(any(seen(replay_depth:)<1).or.any(abs(cursor(replay_depth:)-ending)>tolerance)) return
      identity=replay_interval_history(1);identity%initial=identity%clock_initial
      if(.not.valid_held_mol_record_3d(tree,identity,nv,frame)) return
      if(.not.valid_held_mol_field_3d(tree,identity%clock_internal_energy,1).or. &
        .not.valid_held_mol_field_3d(tree,identity%clock_reaction_source,nv-4)) return
      if(present(heat_release)) then
        if(.not.valid_held_mol_field_3d(tree,identity%clock_heat,1)) return
      endif
      if(replay_depth>0) then
        identity%initial=replay_interval_history(1)%ancestor_start
        if(.not.valid_held_mol_record_3d(tree,identity,nv,frame)) return
        identity%initial=replay_interval_history(1)%ancestor_end
        if(.not.valid_held_mol_record_3d(tree,identity,nv,frame)) return
        if(.not.valid_held_mol_field_3d(tree,identity%ancestor_internal_energy_start,1).or. &
          .not.valid_held_mol_field_3d(tree,identity%ancestor_internal_energy_end,1).or. &
          .not.valid_held_mol_field_3d(tree,identity%ancestor_source_start,nv).or. &
          .not.valid_held_mol_field_3d(tree,identity%ancestor_source_end,nv)) return
        if(.not.allocated(identity%open_ancestor_registers).or. &
          .not.allocated(identity%open_ancestor_held_registers)) return
        if(size(identity%open_ancestor_registers)/=size(tree%nodes).or. &
          size(identity%open_ancestor_held_registers)/=size(tree%nodes)) return
        do c=1,replay_depth
          parent_dt=dt/(identity%fractions_end(c)-identity%fractions_start(c))
          parent_start=replay_time-identity%fractions_start(c)*parent_dt
          if(.not.ieee_is_finite(parent_dt).or..not.ieee_is_finite(parent_start).or.parent_dt<=0) return
          do r=1,size(replay_interval_history)
            expected_start=(replay_interval_history(r)%time-parent_start)/parent_dt
            expected_end=(replay_interval_history(r)%time+replay_interval_history(r)%interval-parent_start)/parent_dt
            if(abs(expected_start-replay_interval_history(r)%fractions_start(c))>32*epsilon(1._dp).or. &
              abs(expected_end-replay_interval_history(r)%fractions_end(c))>32*epsilon(1._dp)) return
          enddo
        enddo
        do p=1,size(tree%nodes)
          d=tree%topology%nodes(p)%shape
          call initialize_state_redist_tree_register_3d(tree%topology,p,nv,6,expected_register,valid, &
            boundary_faces,scale=.5_dp)
          if(.not.valid) return
          call accumulate_state_redist_register_3d(expected_register,identity%open_ancestor_registers(p),nv,d+12,valid)
          if(.not.valid) return
          call initialize_state_redist_tree_register_3d(tree%topology,p,1,6,expected_register,valid, &
            boundary_faces,scale=.5_dp)
          if(.not.valid) return
          call accumulate_state_redist_register_3d(expected_register,identity%open_ancestor_held_registers(p),1,d+12,valid)
          if(.not.valid) return
          if(allocated(identity%open_ancestor_registers(p)%coarse_delta)) then
            if(any(.not.ieee_is_finite(identity%open_ancestor_registers(p)%coarse_delta))) return
          endif
          if(allocated(identity%open_ancestor_registers(p)%fine_delta)) then
            if(any(.not.ieee_is_finite(identity%open_ancestor_registers(p)%fine_delta))) return
          endif
          if(allocated(identity%open_ancestor_held_registers(p)%coarse_delta)) then
            if(any(.not.ieee_is_finite(identity%open_ancestor_held_registers(p)%coarse_delta))) return
          endif
          if(allocated(identity%open_ancestor_held_registers(p)%fine_delta)) then
            if(any(.not.ieee_is_finite(identity%open_ancestor_held_registers(p)%fine_delta))) return
          endif
        enddo
      endif
      good=.true.
    end subroutine validate_prepared_clock_history

    subroutine frozen_auxiliary_source(state,interval,stage,depth,rate,good)
      type(eb_reactive_tree_3d), intent(in) :: state
      real(dp), intent(in) :: interval
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      logical, intent(out) :: good
      good=.false.
    end subroutine frozen_auxiliary_source

    subroutine source_adapter(state,interval,stage,depth,rate,ok)
      type(eb_reactive_tree_3d), intent(in) :: state
      real(dp), intent(in) :: interval
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      logical, intent(out) :: ok
      type(amr_tree_state_3d) :: off
      integer :: p,x,y,z,d(3)
      logical :: valid,agreed
      call auxiliary_source(state,interval,stage,depth,source_time,rate,ok)
      if(.not.routing) return
      if(.not.ok) return
      ok=.false.
      if(.not.allocated(rate%fields)) return
      if(size(rate%fields)/=size(tree%nodes).or.lbound(rate%fields,1)/=1) return
      off=rate;valid=.true.
      do p=1,size(tree%nodes)
        if(allocated(off%fields(p)%values)) off%fields(p)%values=0
      end do
      do p=1,size(tree%nodes)
        if(tree%topology%nodes(p)%level==depth) cycle
        d=tree%topology%nodes(p)%shape
        if(.not.allocated(rate%fields(p)%values)) return
        if(any(shape(rate%fields(p)%values)/=[naux,d]).or.any(lbound(rate%fields(p)%values)/=1)) return
        off%fields(p)%values=0
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(state%nodes(p)%geometry%volume_fraction(x,y,z)<=0.or.covered(tree%topology,p,[x,y,z])) cycle
          off%fields(p)%values(:,x,y,z)=rate%fields(p)%values(:,x,y,z)
        end do
        end do
        end do
        valid=valid.and.all(ieee_is_finite(off%fields(p)%values))
      end do
      if(.not.valid) return
      if(stage==0) then
        off_old=off
      else
        off_new=off
      end if
      ok=.true.

    end subroutine source_adapter
    subroutine build_held_parent_trial(first,weight,result,good)
      integer, intent(in) :: first
      real(dp), intent(in) :: weight
      type(held_replay_result), intent(out) :: result
      logical, intent(out) :: good
      logical :: refused
      real(dp) :: trial_boundary(nv),trial_theta
      refused=.false.;trial_boundary=0;trial_theta=1
      result%tree=candidate;result%stored=held_candidate;result%reaction_source=ir
      result%records=held_history(first:);result%prefix=correction_work
      result%chemical=chemical
      if(present(heat_release)) allocate(result%heat,source=heat)
      if(prepared_redist) then
        allocate(result%wall(nv),result%source(nv))
        if(has_auxiliary) allocate(result%auxiliary(naux))
        call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,held_history(first)%interval, &
          iterations,chemistry,transport_flags,force,rtol,atol,result%tree,extended_geometry,boundary_faces, &
          trial_boundary,result%wall,result%chemical,result%source,trial_theta,good,result%reaction_source, &
          difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature, &
          redistribution=method,target_volume_fraction=redist_target,evaluate_reaction=evaluate_reaction, &
          mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,embedded_wall_temperature=embedded_wall_temperature, &
          embedded_noslip=embedded_noslip,frame=frame,state_redist_max_order=redist_order, &
          evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
        spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=result%heat,integrated_flux=result%flux, &
          internal_energy=result%stored,internal_energy_flux=result%scalar_flux,energy_policy=energy_policy, &
          energy_correction=result%prefix,held_interval_history=result%records, &
          replay_interval_history=held_history(first:),replay_spatial_weight=weight, &
          replay_parent_registers=result%parent_registers,replay_parent_held_registers=result%parent_held_registers, &
          provider_refused=refused,auxiliary_change=result%auxiliary,route_auxiliary=routing, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        observe_mol=observe_mol,observe_mol_sync=observe_mol_sync,observe_mol_reaction=observe_mol_reaction)
      else if(allocated(funding_context_work)) then
        ! Replay the actual recursive clock and pre-payment anchor. A legacy
        ! face-only group cannot migrate pending receipts or open ancestors.
        allocate(result%funding,source=funding_context_work)
        allocate(result%wall(nv),result%source(nv))
        if(has_auxiliary) allocate(result%auxiliary(naux))
        result%wall=0;result%source=0
        if(has_auxiliary) result%auxiliary=0
        call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,held_history(first)%interval, &
          iterations,chemistry,transport_flags,force,rtol,atol,result%tree,extended_geometry,boundary_faces, &
          trial_boundary,result%wall,result%chemical,result%source,trial_theta,good,result%reaction_source, &
          difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature, &
          redistribution=method,target_volume_fraction=redist_target,evaluate_reaction=evaluate_reaction, &
          mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs,embedded_wall_temperature=embedded_wall_temperature, &
          embedded_noslip=embedded_noslip,frame=frame,state_redist_max_order=redist_order, &
          evaluate_spatial_diffusion=evaluate_spatial_diffusion,spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
          constant_coefficients=constant_coefficients,heat_release=result%heat,integrated_flux=result%flux, &
          internal_energy=result%stored,internal_energy_flux=result%scalar_flux,energy_policy=energy_policy, &
          energy_correction=result%prefix,held_interval_history=result%records, &
          native_replay_history=held_history(first:),native_spatial_weight=weight, &
          source_funding_context=result%funding,provider_refused=refused, &
          auxiliary_source=auxiliary_source,auxiliary_change=result%auxiliary,route_auxiliary=routing, &
          spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        observe_mol=observe_mol,observe_mol_sync=observe_mol_sync,observe_mol_reaction=observe_mol_reaction)
      else
      call replay_reactive_eb_mol_held_group_3d(species,reactions,result%tree,result%stored, &
        result%reaction_source,result%prefix,result%records,weight,result%flux,result%scalar_flux, &
        result%chemical,rtol,atol,good,boundary_faces,evaluate_reaction,frame,result%heat,refused)
      endif
      held_replay_provider_refused=held_replay_provider_refused.or.refused
    end subroutine build_held_parent_trial

    subroutine repair_held_parent(first,flux,scalar_flux,good)
      integer, intent(in) :: first
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:),scalar_flux(:)
      logical, intent(out) :: good
      type(held_replay_result) :: best,trial
      type(eb_integrated_flux_3d), allocatable :: new_registration(:)
      real(dp) :: lower,upper,weight,new_chemical(nv),new_wall(nv),new_source(nv),new_auxiliary(naux)
      integer :: attempt,r,p,depth
      logical :: accepted
      good=.false.
      call build_held_parent_trial(first,0._dp,best,accepted)
      if(.not.accepted.or.held_replay_provider_refused) return
      lower=0;upper=1
      call build_held_parent_trial(first,1._dp,trial,accepted)
      if(held_replay_provider_refused) return
      if(accepted) then
        best=trial;lower=1
      else
        do attempt=1,12
          weight=.5_dp*(lower+upper)
          call build_held_parent_trial(first,weight,trial,accepted)
          if(held_replay_provider_refused) return
          if(accepted) then
            best=trial;lower=weight
          else
            upper=weight
          endif
        enddo
      endif
      new_chemical=0;new_wall=0;new_source=0;new_auxiliary=0
      if(present(registered_flux)) call allocate_held_mol_faces_3d(candidate,nv,new_registration)
      do r=1,size(held_history)
        if(r>=first) then
          new_chemical=new_chemical+best%records(r-first+1)%chemical_change
          if(prepared_redist.or.allocated(best%funding)) then
            new_wall=new_wall+best%records(r-first+1)%wall_change
            new_source=new_source+best%records(r-first+1)%source_change
            new_auxiliary=new_auxiliary+best%records(r-first+1)%auxiliary_change
          endif
        else
          new_chemical=new_chemical+held_history(r)%chemical_change
          if(prepared_redist.or.allocated(best%funding)) then
            new_wall=new_wall+held_history(r)%wall_change
            new_source=new_source+held_history(r)%source_change
            new_auxiliary=new_auxiliary+held_history(r)%auxiliary_change
          endif
        endif
        if(.not.present(registered_flux)) cycle
        do p=1,size(tree%nodes)
          if(tree%topology%nodes(p)%level/=held_history(r)%depth) cycle
          if(r>=first) then
            new_registration(p)%x=new_registration(p)%x+best%records(r-first+1)%spatial_flux(p)%x
            new_registration(p)%y=new_registration(p)%y+best%records(r-first+1)%spatial_flux(p)%y
            new_registration(p)%z=new_registration(p)%z+best%records(r-first+1)%spatial_flux(p)%z
          else
            new_registration(p)%x=new_registration(p)%x+held_history(r)%spatial_flux(p)%x
            new_registration(p)%y=new_registration(p)%y+held_history(r)%spatial_flux(p)%y
            new_registration(p)%z=new_registration(p)%z+held_history(r)%spatial_flux(p)%z
          endif
        enddo
      enddo
      if(any(.not.ieee_is_finite(new_chemical))) return
      if(prepared_redist.or.allocated(best%funding)) then
        if(any(.not.ieee_is_finite(new_wall)).or.any(.not.ieee_is_finite(new_source)).or. &
          any(.not.ieee_is_finite(new_auxiliary))) return
      endif
      candidate=best%tree;held_candidate=best%stored;ir=best%reaction_source
      correction_work=best%prefix;chemical=new_chemical;held_history(first:)=best%records
      held_history(first:)%replayed=.true.
      if(allocated(best%funding)) then
        call move_eb_held_funding_context_3d(best%funding,funding_context_work)
        wall=new_wall;external=new_source;auxiliary=new_auxiliary
      endif
      if(present(heat_release)) heat=best%heat
      if(present(registered_flux)) call move_alloc(new_registration,registration)
      depth=held_history(first)%depth;deferred_record(depth:)=0
      if(prepared_redist) then
        wall=new_wall;external=new_source;auxiliary=new_auxiliary
        if(depth>0) then
          interface_registers(depth-1)%fields=best%parent_registers
          held_interface_registers(depth-1)%fields=best%parent_held_registers
        endif
        do r=depth,levels-2
          do p=1,size(tree%nodes)
            if(allocated(interface_registers(r)%fields(p)%coarse_delta)) interface_registers(r)%fields(p)%coarse_delta=0
            if(allocated(interface_registers(r)%fields(p)%fine_delta)) interface_registers(r)%fields(p)%fine_delta=0
            if(allocated(held_interface_registers(r)%fields(p)%coarse_delta)) &
              held_interface_registers(r)%fields(p)%coarse_delta=0
            if(allocated(held_interface_registers(r)%fields(p)%fine_delta)) &
              held_interface_registers(r)%fields(p)%fine_delta=0
          enddo
        enddo
      endif
      do p=1,size(tree%nodes)
        if(tree%topology%nodes(p)%level<depth) cycle
        if(prepared_redist.or.allocated(best%funding)) then
          do r=first,size(held_history)
            if(tree%topology%nodes(p)%level/=held_history(r)%depth) cycle
            starts%nodes(p)=held_history(r)%initial%nodes(p);ends%nodes(p)=held_history(r)%local_tree%nodes(p)
            held_starts%fields(p)=held_history(r)%initial_internal_energy%fields(p)
            held_ends%fields(p)=held_history(r)%local_internal_energy%fields(p)
            if(allocated(best%funding)) then
              source0%fields(p)=held_history(r)%old_predictor%fields(p)
              source1%fields(p)=held_history(r)%new_predictor%fields(p)
            endif
          enddo
        else
          ends%nodes(p)=candidate%nodes(p);held_ends%fields(p)=held_candidate%fields(p)
        endif
      enddo
      call move_alloc(best%flux,flux);call move_alloc(best%scalar_flux,scalar_flux)
      theta=min(theta,lower);good=.true.
    end subroutine repair_held_parent

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
      good=local_ok;if(.not.good) return
        endif
        if(allocated(off_new%fields)) then
          call copy_godunov_field_3d(off_new,saved%routed_auxiliary_new_rate,local_ok)
      good=local_ok;if(.not.good) return
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
        good=local_ok;if(.not.good) return
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
      logical :: admitted
      good=.false.
      if(depth==levels-1) then
        good=.true.;return
      endif
      pressure=held_candidate
      do p=1,size(tree%nodes)
        pressure%fields(p)%values=0
        if(tree%topology%nodes(p)%level<depth) cycle
        m=tree%topology%nodes(p)%shape
        do z=1,m(3)
        do y=1,m(2)
        do x=1,m(1)
          if(candidate%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
          call reactive_conserved_to_primitive(species,candidate%nodes(p)%state(:,x,y,z), &
            candidate%nodes(p)%temperature(x,y,z),q,t,cs,admitted, &
            internal_energy_density=held_candidate%fields(p)%values(1,x,y,z))
          if(.not.admitted) return
          pressure%fields(p)%values(1,x,y,z)=q(5)
        enddo
        enddo
        enddo
      enddo
      call average_down_reactive_eb_tree_field_3d(candidate,pressure,admitted,minimum_depth=depth)
      if(.not.admitted) return
      do p=1,size(tree%nodes)
        if(tree%topology%nodes(p)%level/=depth) cycle
        m=tree%topology%nodes(p)%shape
        do z=1,m(3)
        do y=1,m(2)
        do x=1,m(1)
          if(candidate%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
          if(.not.covered(tree%topology,p,[x,y,z])) cycle
          call reactive_conserved_to_primitive(species,candidate%nodes(p)%state(:,x,y,z), &
            candidate%nodes(p)%temperature(x,y,z),q,t,cs,admitted, &
            internal_energy_density=held_candidate%fields(p)%values(1,x,y,z))
          if(.not.admitted) return
          target=pressure%fields(p)%values(1,x,y,z)
          ! Re-entering a captured clock must not accumulate roundoff resets.
          if(abs(q(5)-target)<=64*epsilon(1._dp)*max(abs(q(5)),abs(target))) cycle
          call reconstruct_eb_pressure_state_3d(species,candidate%nodes(p)%state(:,x,y,z), &
            held_candidate%fields(p)%values(1,x,y,z),candidate%nodes(p)%temperature(x,y,z),target,delta,admitted)
          if(.not.admitted) return
        enddo
        enddo
        enddo
      enddo
      good=.true.
    end subroutine prepare_clock_pressure

    recursive subroutine advance_level(level,time,interval,flux,physical_flux,scalar_flux,good)
      type(amr_tree_state_3d), allocatable :: diagnostic_raw,diagnostic_applied
      integer, intent(in) :: level
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
      type(eb_integrated_flux_3d), allocatable :: sync_flux(:),paired_flux(:)
      type(eb_integrated_flux_3d), allocatable :: registered_replacement(:)
      type(eb_integrated_flux_3d), allocatable :: step_registration(:)
      type(amr_tree_state_3d) :: op,np,correction,redistributed,registered_increment
      type(amr_tree_state_3d) :: off_transfer
      type(eb_reactive_tree_3d) :: off_ledger
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

      type(amr_tree_state_3d), allocatable :: transfer
      type(state_redist_register_3d), allocatable :: stage_registers(:),held_stage_registers(:)
      type(state_redist_register_3d) :: part
      type(eb_reactive_tree_3d) :: before_sync
      type(eb_mol_held_interval_3d), allocatable :: step_history,step_replay
      procedure(eb_sdc_auxiliary_source), pointer :: prepared_auxiliary
      type(eb_mol_held_interval_3d) :: clock_record
      integer :: history_first
      real(dp) :: fstart(level),fend(level),b(nv),w(nv),ch(nv),ext(nv),limiter,elapsed,sub_dt
      integer :: p,d(3),step,target,bucket
      integer :: profile_status
      real(dp) :: clock_start_reset,clock_end_reset,scalar_reset_change
      logical :: accepted,old_only,paired_closed,replayed_closure,step_refused
      real(dp) :: clock_tolerance
      good=.false.;paired_closed=.false.;replayed_closure=.false.;step_refused=.false.
      if(selected_pressure) then
        call prepare_clock_pressure(level,accepted)
        if(.not.accepted) return
      endif
      partial_source_replay=.false.;complete_routed_replay=.false.
      clock_start_reset=0;clock_end_reset=0;scalar_reset_change=0
      if(replaying_prepared) then
        if(replay_cursor>=size(replay_interval_history)) return
        replay_cursor=replay_cursor+1
        clock_tolerance=32*epsilon(1._dp)*max(abs(time+interval),tiny(1._dp))
        if(replay_interval_history(replay_cursor)%depth/=level.or. &
          abs(replay_interval_history(replay_cursor)%time-time)>clock_tolerance.or. &
          abs(replay_interval_history(replay_cursor)%interval-interval)>clock_tolerance) return
        allocate(step_replay,source=replay_interval_history(replay_cursor))
        do p=1,size(tree%nodes)
          if(tree%topology%nodes(p)%level==level) step_replay%initial_reaction_source%fields(p)=ir%fields(p)
        enddo
      endif
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
        if(.not.(saved%depth==level.and.saved%time==time.and.saved%interval==interval)) return
        end associate
        allocate(step_replay,source=native_replay_history(native_cursor))
        do p=1,size(tree%nodes)
          if(tree%topology%nodes(p)%level==level) step_replay%initial_reaction_source%fields(p)=ir%fields(p)
        enddo
      endif
      if(held_replay_enabled) deferred_record(level)=0
      old_only=.false.
      if(present(mol_mode)) old_only=mol_mode.and..not.chemistry.and.iterations>1
      if(present(godunov_replay_history)) then
        if(godunov_replay_cursor>=godunov_replay_count) return
        godunov_replay_cursor=godunov_replay_cursor+1
        if(godunov_replay_cursor<=size(godunov_replay_history)) then
        allocate(step_replay_godunov,source=godunov_replay_history(godunov_replay_cursor),stat=godunov_status)
        if(godunov_status/=0) return
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
          call copy_godunov_field_3d(step_replay_godunov%routed_auxiliary_old_rate,off_old,accepted)
          if(.not.accepted) return
        endif
        if(routing.and.allocated(step_replay_godunov%routed_auxiliary_new_rate%fields)) then
          call copy_godunov_field_3d(step_replay_godunov%routed_auxiliary_new_rate,off_new,accepted)
          if(.not.accepted) return
        endif
      endif
      if(counts(level+1)==huge(1_int64)) return
      if(godunov_capture) then
        if(godunov_cursor>=size(godunov_records)) return
        godunov_cursor=godunov_cursor+1;godunov_index=godunov_cursor
        godunov_active(level)=godunov_index+godunov_replay_base
        allocate(step_godunov,step_source_recovery,source_prefix_ready,stat=godunov_status)
        if(godunov_status/=0) return
        godunov_clock%parent_record=0
        if(level>0) godunov_clock%parent_record=godunov_active(level-1)
        call capture_godunov_clock(level,godunov_clock,accepted)
        if(.not.accepted) return
      endif
      if(allocated(funding_context_work)) then
        call copy_eb_held_funding_context_3d(tree%topology,funding_geometry,nv,funding_local,funding_context_work, &
          funding_clock_context,accepted)
        if(.not.accepted) return
        if(held_capture_enabled) then
          clock_record%clock_initial=candidate;clock_record%clock_internal_energy=held_candidate
          clock_record%clock_reaction_source=ir;clock_record%clock_energy_correction=correction_work
          if(present(heat_release)) clock_record%clock_heat=heat
          clock_record%clock_source_funding_context=funding_clock_context
        endif
        call take_eb_held_funding_target_3d(tree%topology,funding_geometry,nv,funding_local,level, &
          funding_context_work%pending,funding_payment,funding_queue_next,accepted)
        if(.not.accepted) return
        call apply_native_funding_receipt(funding_payment,accepted);if(.not.accepted) return
        call move_alloc(funding_queue_next,funding_context_work%pending)
      endif
      if(allocated(correction_work)) then
        call reset_held_eb_tree_family_3d(species,candidate,held_candidate,reset_difference,held_policy,level, &
          reset_change,accepted,frame, &
          internal_energy_change=clock_start_reset,reset_level=merge(-2,-1,godunov_mol.and.chemistry.and.method=='none'))
        if(.not.accepted) return
        correction_work=correction_work+reset_change
      endif
      if(godunov_capture) then
        godunov_clock%held_clock_start_change=clock_start_reset
        if(godunov_index==1) godunov_clock%held_clock_start_change=clock_start_reset+initial_clock_reset
      endif
      clock_start(level+1)=time;clock_dt(level+1)=interval
      if(held_capture_enabled) then
        allocate(clock_record%clock_start(levels),clock_record%clock_interval(levels),stat=profile_status)
        if(profile_status/=0) return
        clock_record%clock_start=clock_start;clock_record%clock_interval=clock_dt
      endif
      if(replaying_prepared) then
        fstart=0;fend=0
        fstart(:replay_depth)=step_replay%fractions_start(:replay_depth)
        fend(:replay_depth)=step_replay%fractions_end(:replay_depth)
        fstart(replay_depth+1:level)=max(0._dp,min(1._dp, &
          (time-clock_start(replay_depth+1:level))/clock_dt(replay_depth+1:level)))
        fend(replay_depth+1:level)=max(fstart(replay_depth+1:level),min(1._dp, &
          (time+interval-clock_start(replay_depth+1:level))/clock_dt(replay_depth+1:level)))
      else
      fstart=max(0._dp,min(1._dp,(time-clock_start(:level))/clock_dt(:level)))
      fend=max(fstart,min(1._dp,(time+interval-clock_start(:level))/clock_dt(:level)))
      endif
      do p=1,size(tree%nodes)
        if(tree%topology%nodes(p)%level==level) starts%nodes(p)=candidate%nodes(p)
        if(present(internal_energy)) then
          if(tree%topology%nodes(p)%level==level) held_starts%fields(p)=held_candidate%fields(p)
        endif
      end do
      if(held_capture_enabled) then
        allocate(step_history)
        clock_record%held_clock_start_change=clock_start_reset
        if(size(held_history)==0) then
          clock_record%held_clock_start_change=clock_start_reset+initial_clock_reset
          if(replaying_prepared) then
            call compute_reactive_eb_held_projection(0._dp,replay_interval_history(1)%held_clock_start_change,0._dp, &
              clock_record%held_clock_start_change,scalar_reset_change,accepted)
            if(.not.accepted) return
            clock_record%held_clock_start_change=scalar_reset_change
          endif
        endif
        if(replaying_native.and.size(held_history)==0) then
          call compute_reactive_eb_held_projection(0._dp,native_replay_history(1)%held_clock_start_change,0._dp, &
            clock_record%held_clock_start_change,scalar_reset_change,accepted)
          if(.not.accepted) return
          clock_record%held_clock_start_change=scalar_reset_change
        endif
        if(.not.allocated(funding_context_work)) then
        clock_record%clock_initial=candidate;clock_record%clock_internal_energy=held_candidate
        clock_record%clock_reaction_source=ir;clock_record%clock_energy_correction=correction_work
        if(present(heat_release)) clock_record%clock_heat=heat
        endif
        if(prepared_redist.and.level>0) then
          clock_record%open_ancestor_registers=interface_registers(level-1)%fields
          clock_record%open_ancestor_held_registers=held_interface_registers(level-1)%fields
        endif
      endif
      if(method/='none'.and..not.prepared_redist) allocate(transfer)
      nullify(prepared_auxiliary)
      if(has_auxiliary) then
        allocate(step_aux(naux))
        if(replaying_prepared.or.step_replaying_native) then
          prepared_auxiliary=>frozen_auxiliary_source
        else
          prepared_auxiliary=>source_adapter
        endif
      endif
      if(present(internal_energy_pressure_work_change).or.godunov_capture) &
        allocate(step_pressure_work,source=0._dp)
      if(allocated(funding_context_work)) then
        call read_eb_held_funding_open_faces_3d(tree%topology,funding_geometry,nv,funding_local,level, &
          funding_context_work,funding_faces,funding_scalar_faces,accepted)
        if(.not.accepted) return
        call assemble_reactive_eb_held_funding_debt_3d(species,candidate,held_candidate,funding_faces, &
          funding_scalar_faces,level,funding_debt,accepted,boundary_faces,frame)
        if(.not.accepted) return
        allocate(funding_stock,funding_receipt,stat=funding_status);if(funding_status/=0) return
        call read_eb_held_funding_stock_3d(tree%topology,funding_geometry,nv,funding_local,level, &
          funding_context_work%pending,funding_stock,accepted,open_debt=funding_debt)
        if(.not.accepted) return
      endif
      source_time=time
      native_retry_enabled=(held_replay_enabled.or.(native_prefix.and..not.step_replaying_native)).and. &
        present(energy_correction).and.method=='none'.and..not.routing.and. &
        .not.present(evaluate_spatial_diffusion)
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
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,interval,iterations,chemistry, &
        transport_flags,force,rtol,atol,candidate,extended_geometry,boundary_faces,b,w,ch,ext,limiter,accepted,ir, &
        difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature, &
        level,starts,ends,fstart,fend,source0,source1,flux,op,np,method, &
        evaluate_reaction=evaluate_reaction,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs, &
        frame=frame,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
        spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat, &
        embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip, &
        registered_flux=step_registration,internal_energy=held_candidate,internal_energy_flux=step_scalar_flux, &
        ancestor_internal_energy_start=held_starts,ancestor_internal_energy_end=held_ends, &
        target_volume_fraction=redist_target,state_redist_max_order=redist_order, &
        state_redist_registers=stage_registers,held_state_redist_registers=held_stage_registers, &
        energy_policy=energy_policy,energy_correction=correction_work,held_interval_history=step_history, &
        replay_interval=step_replay,prepared_spatial_weight=replay_weight, &
        auxiliary_source=prepared_auxiliary,auxiliary_change=step_aux, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        internal_energy_pressure_work_change=step_pressure_work, &
        godunov_interval_history=step_godunov,godunov_replay_interval=step_replay_godunov, &
        godunov_replay_partial=partial_source_replay,godunov_source_recovery=step_source_recovery, &
        godunov_source_recovery_ready=source_prefix_ready, &
        spatial_flux_weight=spatial_flux_weight,provider_refused=step_refused, &
        source_funding_stock=funding_stock,source_funding_impulse=funding_receipt, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        observe_mol=observe_mol,diagnostic_time=time,observe_mol_reaction=observe_mol_reaction)
      else
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,interval,iterations,chemistry, &
        transport_flags,force,rtol,atol,candidate,extended_geometry,boundary_faces,b,w,ch,ext,limiter,accepted,ir, &
        difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature, &
        level,starts,ends,fstart,fend,source0,source1,flux,op,np,method, &
        evaluate_reaction=evaluate_reaction,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs, &
        frame=frame,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
        spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat, &
        embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip, &
        registered_flux=step_registration,internal_energy=held_candidate,internal_energy_flux=step_scalar_flux, &
        ancestor_internal_energy_start=held_starts,ancestor_internal_energy_end=held_ends, &
        target_volume_fraction=redist_target,state_redist_max_order=redist_order, &
        state_redist_registers=stage_registers,held_state_redist_registers=held_stage_registers, &
        energy_policy=energy_policy,energy_correction=correction_work,held_interval_history=step_history, &
        replay_interval=step_replay,prepared_spatial_weight=replay_weight, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        internal_energy_pressure_work_change=step_pressure_work, &
        godunov_interval_history=step_godunov,godunov_replay_interval=step_replay_godunov, &
        godunov_replay_partial=partial_source_replay,godunov_source_recovery=step_source_recovery, &
        godunov_source_recovery_ready=source_prefix_ready, &
        spatial_flux_weight=spatial_flux_weight,provider_refused=step_refused, &
        source_funding_stock=funding_stock,source_funding_impulse=funding_receipt, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        observe_mol=observe_mol,diagnostic_time=time,observe_mol_reaction=observe_mol_reaction)
      endif
      held_replay_provider_refused=held_replay_provider_refused.or.step_refused
      else
      if(has_auxiliary) then
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,interval,iterations,chemistry, &
        transport_flags,force,rtol,atol,candidate,extended_geometry,boundary_faces,b,w,ch,ext,limiter,accepted,ir, &
        difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature, &
        level,starts,ends,fstart,fend,source0,source1,flux,op,np,method, &
        evaluate_reaction=evaluate_reaction,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs, &
        frame=frame,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
        spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat, &
        embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip, &
        registered_flux=step_registration,internal_energy=held_candidate,internal_energy_flux=step_scalar_flux, &
        ancestor_internal_energy_start=held_starts,ancestor_internal_energy_end=held_ends, &
        energy_policy=energy_policy,energy_correction=correction_work,held_interval_history=step_history, &
        auxiliary_source=prepared_auxiliary,auxiliary_change=step_aux, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        internal_energy_pressure_work_change=step_pressure_work, &
        godunov_interval_history=step_godunov,godunov_replay_interval=step_replay_godunov, &
        godunov_replay_partial=partial_source_replay,godunov_source_recovery=step_source_recovery, &
        godunov_source_recovery_ready=source_prefix_ready, &
        spatial_flux_weight=spatial_flux_weight,provider_refused=step_refused, &
        source_funding_stock=funding_stock,source_funding_impulse=funding_receipt, &
        replay_interval=step_replay,freeze_native_external_sources=step_replaying_native, &
        native_spatial_weight=common_native_weight, &
        native_source_prefix=native_prefix_output,native_source_prefix_ready=native_prefix_ready, &
        native_source_replay=native_prefix_input, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        observe_mol=observe_mol,diagnostic_time=time,observe_mol_reaction=observe_mol_reaction)
      else
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,interval,iterations,chemistry, &
        transport_flags,force,rtol,atol,candidate,extended_geometry,boundary_faces,b,w,ch,ext,limiter,accepted,ir, &
        difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature, &
        level,starts,ends,fstart,fend,source0,source1,flux,op,np,method, &
        evaluate_reaction=evaluate_reaction,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs, &
        frame=frame,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
        spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat, &
        embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip, &
        registered_flux=step_registration,internal_energy=held_candidate,internal_energy_flux=step_scalar_flux, &
        ancestor_internal_energy_start=held_starts,ancestor_internal_energy_end=held_ends, &
        energy_policy=energy_policy,energy_correction=correction_work,held_interval_history=step_history, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        internal_energy_pressure_work_change=step_pressure_work, &
        godunov_interval_history=step_godunov,godunov_replay_interval=step_replay_godunov, &
        godunov_replay_partial=partial_source_replay,godunov_source_recovery=step_source_recovery, &
        godunov_source_recovery_ready=source_prefix_ready, &
        spatial_flux_weight=spatial_flux_weight,provider_refused=step_refused, &
        source_funding_stock=funding_stock,source_funding_impulse=funding_receipt, &
        replay_interval=step_replay,freeze_native_external_sources=step_replaying_native, &
        native_spatial_weight=common_native_weight, &
        native_source_prefix=native_prefix_output,native_source_prefix_ready=native_prefix_ready, &
        native_source_replay=native_prefix_input, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        observe_mol=observe_mol,diagnostic_time=time,observe_mol_reaction=observe_mol_reaction)
      endif
      endif
      else if(prepared_redist) then
      if(present(auxiliary_source)) then
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,interval,iterations,chemistry, &
        transport_flags,force,rtol,atol,candidate,extended_geometry,boundary_faces,b,w,ch,ext,limiter,accepted,ir, &
        difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature, &
        level,starts,ends,fstart,fend,source0,source1,flux,op,np,method,transfer,redist_target,source_adapter,step_aux, &
          evaluate_reaction=evaluate_reaction,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs, &
          embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip,frame=frame, &
      state_redist_max_order=redist_order,state_redist_registers=stage_registers, &
          evaluate_spatial_diffusion=evaluate_spatial_diffusion,spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat,registered_flux=step_registration, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        observe_mol=observe_mol,diagnostic_time=time,observe_mol_reaction=observe_mol_reaction)
      else
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,interval,iterations,chemistry, &
        transport_flags,force,rtol,atol,candidate,extended_geometry,boundary_faces,b,w,ch,ext,limiter,accepted,ir, &
        difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature, &
        level,starts,ends,fstart,fend,source0,source1,flux,op,np,method,transfer,redist_target, &
          evaluate_reaction=evaluate_reaction,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs, &
          embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip,frame=frame, &
      state_redist_max_order=redist_order,state_redist_registers=stage_registers, &
          evaluate_spatial_diffusion=evaluate_spatial_diffusion,spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat,registered_flux=step_registration, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        observe_mol=observe_mol,diagnostic_time=time,observe_mol_reaction=observe_mol_reaction)
      end if
      else
      if(has_auxiliary) then
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,interval,iterations,chemistry, &
        transport_flags,force,rtol,atol,candidate,extended_geometry,boundary_faces,b,w,ch,ext,limiter,accepted,ir, &
        difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature, &
        level,starts,ends,fstart,fend,source0,source1,flux,op,np,method,transfer,redist_target,source_adapter,step_aux, &
          evaluate_reaction=evaluate_reaction,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs, &
          embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip,frame=frame, &
      state_redist_max_order=redist_order,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
          spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat,registered_flux=step_registration, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        observe_mol=observe_mol,diagnostic_time=time,observe_mol_reaction=observe_mol_reaction)
      else
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,interval,iterations,chemistry, &
        transport_flags,force,rtol,atol,candidate,extended_geometry,boundary_faces,b,w,ch,ext,limiter,accepted,ir, &
        difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature, &
        level,starts,ends,fstart,fend,source0,source1,flux,op,np,method,transfer,redist_target, &
          evaluate_reaction=evaluate_reaction,mol_mode=mol_mode,mol_lax_friedrichs=mol_lax_friedrichs, &
          embedded_wall_temperature=embedded_wall_temperature,embedded_noslip=embedded_noslip,frame=frame, &
      state_redist_max_order=redist_order,evaluate_spatial_diffusion=evaluate_spatial_diffusion, &
          spatial_diffusion_halo_width=spatial_diffusion_halo_width, &
        constant_coefficients=constant_coefficients,heat_release=heat,registered_flux=step_registration, &
        spatial_diffusion_is_source=external_diffusion_source,auxiliary_has_internal_energy=held_auxiliary, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        observe_mol=observe_mol,diagnostic_time=time,observe_mol_reaction=observe_mol_reaction)
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
      held_replay_provider_refused=held_replay_provider_refused.or.step_refused
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
        if(.not.ieee_is_finite(pressure_work_change)) return
      endif
      if(allocated(funding_context_work)) then
        call queue_eb_held_funding_receipt_3d(tree%topology,funding_geometry,nv,funding_local,level, &
          funding_context_work%pending,funding_receipt,funding_queue_next,accepted)
        if(.not.accepted) return
        call move_alloc(funding_queue_next,funding_context_work%pending)
      endif
      if(held_capture_enabled) then
        step_history%time=time
        step_history%old_predictor=op;step_history%new_predictor=np
        if(replaying_prepared) then
          step_history%old_flux=step_replay%old_flux
          step_history%old_internal_energy_flux=step_replay%old_internal_energy_flux
        endif
        call move_eb_held_funding_context_3d(clock_record%clock_source_funding_context, &
          step_history%clock_source_funding_context)
        step_history%clock_initial=clock_record%clock_initial
        step_history%clock_internal_energy=clock_record%clock_internal_energy
        step_history%clock_reaction_source=clock_record%clock_reaction_source
        step_history%clock_energy_correction=clock_record%clock_energy_correction
        call move_alloc(clock_record%clock_start,step_history%clock_start)
        call move_alloc(clock_record%clock_interval,step_history%clock_interval)
        call compute_reactive_eb_held_projection(0._dp,clock_record%held_clock_start_change,0._dp, &
          step_history%initial_internal_energy_change,step_history%held_clock_start_change,accepted)
        if(.not.accepted) return
        step_history%clock_heat=clock_record%clock_heat
        call move_alloc(clock_record%open_ancestor_registers,step_history%open_ancestor_registers)
        call move_alloc(clock_record%open_ancestor_held_registers,step_history%open_ancestor_held_registers)
        history_first=size(held_history)+1;held_history=[held_history,step_history]
      endif
      if(present(registered_flux)) then
        if(.not.allocated(registration)) then
          allocate(registration(size(tree%nodes)))
          do p=1,size(registration)
            d=tree%topology%nodes(p)%shape
            allocate(registration(p)%x(nv,0:d(1),d(2),d(3)),registration(p)%y(nv,d(1),0:d(2),d(3)), &
              registration(p)%z(nv,d(1),d(2),0:d(3)),source=0._dp)
          enddo
        endif
        do p=1,size(registration)
          if(tree%topology%nodes(p)%level/=level) then
            if(.not.allocated(step_registration(p)%x).and..not.allocated(step_registration(p)%y).and. &
              .not.allocated(step_registration(p)%z)) cycle
          endif
          if(.not.valid_eb_integrated_flux_3d(step_registration(p),candidate%nodes(p)%geometry,nv)) return
          registration(p)%x=registration(p)%x+step_registration(p)%x
          registration(p)%y=registration(p)%y+step_registration(p)%y
          registration(p)%z=registration(p)%z+step_registration(p)%z
          if(.not.all(ieee_is_finite(registration(p)%x)).or..not.all(ieee_is_finite(registration(p)%y)).or. &
            .not.all(ieee_is_finite(registration(p)%z))) return
        enddo
      endif
      if(prepared_redist) then
        do p=1,size(tree%nodes)
          if(tree%topology%nodes(p)%level/=level) cycle
          d=tree%topology%nodes(p)%shape
          if(level+1<levels) then
            part=stage_registers(p)
            if(allocated(part%fine_delta)) part%fine_delta=0
            call accumulate_state_redist_register_3d(interface_registers(level)%fields(p),part,nv,d+12,accepted)
            if(.not.accepted) return
            if(present(internal_energy)) then
              part=held_stage_registers(p)
              if(allocated(part%fine_delta)) part%fine_delta=0
              call accumulate_state_redist_register_3d(held_interface_registers(level)%fields(p),part,1,d+12,accepted)
              if(.not.accepted) return
            endif
          endif
          if(level>0) then
            part=stage_registers(p)
            if(allocated(part%coarse_delta)) part%coarse_delta=0
            call accumulate_state_redist_register_3d(interface_registers(level-1)%fields(p),part,nv,d+12,accepted)
            if(.not.accepted) return
            if(present(internal_energy)) then
              part=held_stage_registers(p)
              if(allocated(part%coarse_delta)) part%coarse_delta=0
              call accumulate_state_redist_register_3d(held_interface_registers(level-1)%fields(p),part,1,d+12,accepted)
              if(.not.accepted) return
            endif
          endif
        enddo
      else if(method/='none') then
        do p=1,size(tree%nodes)
          target=tree%topology%nodes(p)%level;bucket=min(level,target)
          pending(bucket)%fields(p)%values=pending(bucket)%fields(p)%values+transfer%fields(p)%values
        end do
      end if
      if(routing) then
        if(complete_routed_replay) then
          call copy_godunov_field_3d(step_replay_godunov%routed_auxiliary_impulse,off_transfer,accepted)
          if(.not.accepted) return
        else if(replaying_prepared) then
          off_transfer=step_replay%routed_source
        else
        off_transfer=off_old
        do p=1,size(tree%nodes)
          target=tree%topology%nodes(p)%level;bucket=min(level,target)
          if(old_auxiliary_only) then
            off_transfer%fields(p)%values=.5_dp*interval*off_old%fields(p)%values
          else
            off_transfer%fields(p)%values=.5_dp*interval*(off_old%fields(p)%values+off_new%fields(p)%values)
          endif
        enddo
        endif
        do p=1,size(tree%nodes)
          target=tree%topology%nodes(p)%level;bucket=min(level,target)
          auxiliary_pending(bucket)%fields(p)%values=auxiliary_pending(bucket)%fields(p)%values+ &
            off_transfer%fields(p)%values
        enddo
        if(godunov_capture) then
          call copy_godunov_field_3d(off_transfer,godunov_records(godunov_index)%routed_auxiliary_impulse,accepted)
          if(.not.accepted) return
        endif
        if(held_capture_enabled) then
          held_history(history_first)%routed_auxiliary=.true.
          held_history(history_first)%routed_source=off_transfer
        endif
        off_ledger=candidate
        do p=1,size(tree%nodes)
          off_ledger%nodes(p)%state=off_transfer%fields(p)%values
        end do
        call integrate_reactive_eb_tree_3d(off_ledger,off_inventory,accepted)
        if(.not.accepted) return
        step_aux=step_aux+off_inventory
        if(godunov_capture) godunov_records(godunov_index)%auxiliary_change=step_aux
        if(held_capture_enabled) held_history(history_first)%auxiliary_change=step_aux
      end if
      if(has_auxiliary) then
        auxiliary=auxiliary+step_aux
        if(.not.all(ieee_is_finite(auxiliary))) return
      end if
      counts(level+1)=counts(level+1)+1
      wall=wall+w;chemical=chemical+ch;external=external+ext;theta=min(theta,limiter)
      do p=1,size(tree%nodes)
        if(tree%topology%nodes(p)%level==level) then
          ends%nodes(p)=candidate%nodes(p)
          if(present(internal_energy)) held_ends%fields(p)=held_candidate%fields(p)
          source0%fields(p)%values=op%fields(p)%values;source1%fields(p)%values=np%fields(p)%values
        else
          ! An inactive native clock may carry a conservative donor face.
          ! Initialize only an absent operator; retain an accepted loan family.
          call complete_interval_faces(flux(p),candidate%nodes(p)%geometry,nv,accepted)
          if(.not.accepted) return
          if(present(internal_energy)) then
            call complete_interval_faces(step_scalar_flux(p),candidate%nodes(p)%geometry,1,accepted)
            if(.not.accepted) return
          endif
        end if
      end do
      physical_flux=flux
      if(present(internal_energy)) scalar_flux=step_scalar_flux
      if(present(mol_mode)) then
        if(mol_mode.and..not.chemistry.and.iterations>1) then
          do p=1,size(flux)
            if(tree%topology%nodes(p)%level/=level) cycle
            flux(p)=step_registration(p)
          enddo
        endif
      endif
      if(allocated(funding_context_work)) then
        call set_eb_held_funding_clock_faces_3d(tree%topology,funding_geometry,nv,funding_local,level, &
          physical_flux,scalar_flux,funding_context_work,funding_next,accepted)
        if(.not.accepted) return
        call move_eb_held_funding_context_3d(funding_next,funding_context_work)
      endif
      if(level+1<levels) then
        elapsed=0
        do step=1,ratios(level+2)
          sub_dt=interval/real(ratios(level+2),dp)
          if(step==ratios(level+2)) sub_dt=interval-elapsed
          if(sub_dt<=0.or.time+elapsed+sub_dt<=time+elapsed) return
          call advance_level(level+1,time+elapsed,sub_dt,fine_flux,fine_physical_flux,fine_scalar_flux,accepted)
          if(.not.accepted) then
            if(godunov_automatic_repair) go to 880
            if(.not.held_replay_enabled) return
            if(held_replay_provider_refused.or.step/=ratios(level+2).or.deferred_record(level+1)==0) return
            call repair_held_parent(history_first,flux,scalar_flux,accepted)
            if(.not.accepted) then
              if(.not.held_replay_provider_refused) deferred_record(level)=history_first
              return
            endif
            physical_flux=flux;replayed_closure=.true.
            exit
          endif
          do p=1,size(tree%nodes)
            ! Descendant clocks can also fund an exposed ancestor donor.
            ! Its returned faces are increments to this open native family.
            if(tree%topology%nodes(p)%level<=level) then
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
        call take_eb_held_funding_bucket_3d(tree%topology,funding_geometry,nv,funding_local,level, &
          funding_context_work%pending,funding_payment,funding_queue_next,accepted)
        if(.not.accepted) return
        call apply_native_funding_receipt(funding_payment,accepted);if(.not.accepted) return
        call move_alloc(funding_queue_next,funding_context_work%pending)
      endif
      if(prepared_redist) then
        if(replayed_closure) then
          accepted=.true.;paired_closed=.true.
        else if(present(internal_energy).and.allocated(correction_work)) then
          if(old_only) then
            paired_flux=physical_flux
          else
            paired_flux=flux
          endif
          if(level+1<levels) then
            call close_held_prepared_eb_tree_3d(species,candidate,held_candidate,paired_flux,scalar_flux,held_policy, &
              correction_work,level,accepted,boundary_faces,frame,register_geometry, &
              interface_registers(level)%fields,held_interface_registers(level)%fields, &
          internal_energy_change=clock_end_reset)
          else
            call close_held_prepared_eb_tree_3d(species,candidate,held_candidate,paired_flux,scalar_flux,held_policy, &
              correction_work,level,accepted,boundary_faces,frame, &
          extended_geometry=extended_geometry,internal_energy_change=clock_end_reset)
          endif
          if(.not.accepted.and.godunov_automatic_repair) go to 880
          if(.not.accepted.and.held_replay_enabled) then
            call repair_held_parent(history_first,flux,scalar_flux,accepted)
            if(.not.accepted.and.godunov_automatic_repair) go to 880
            if(.not.accepted) then
              if(.not.held_replay_provider_refused) deferred_record(level)=history_first
              return
            endif
            physical_flux=flux;replayed_closure=.true.
          endif
          paired_closed=accepted
          if(accepted.and..not.replayed_closure) then
            if(old_only) then
              call move_alloc(paired_flux,physical_flux)
              call assemble_reactive_eb_tree_reflux_3d(species,candidate,flux,registered_increment, &
                registered_replacement,accepted,level,min(level+1,levels-1),boundary_faces,frame)
              if(.not.accepted.and.godunov_automatic_repair) go to 880
              if(.not.accepted) return
              call move_alloc(registered_replacement,flux)
            else
              call move_alloc(paired_flux,flux)
            endif
          endif
        else if(old_only) then
          ! The pinned MOL registration can contain only the old half while
          ! the accepted gas state contains both physical half steps. Reflux
          ! the accepted state with its physical transfer and retain a
          ! separately reconciled registered transfer for the parent ledger.
          if(level+1<levels) then
            call synchronize_reactive_eb_tree_register_flux_3d(species,candidate,physical_flux,register_geometry, &
              interface_registers(level)%fields,level,accepted,boundary_faces,frame)
          else
            call synchronize_reactive_eb_tree_flux_3d(species,candidate,physical_flux,accepted,level,level, &
              boundary_faces,frame)
          endif
          if(.not.accepted.and.godunov_automatic_repair) go to 880
          if(.not.accepted) return
          call assemble_reactive_eb_tree_reflux_3d(species,candidate,flux,registered_increment, &
            registered_replacement,accepted,level,min(level+1,levels-1),boundary_faces,frame)
          if(.not.accepted.and.godunov_automatic_repair) go to 880
          if(.not.accepted) return
          call move_alloc(registered_replacement,flux)
        else
          if(level+1<levels) then
            call synchronize_reactive_eb_tree_register_flux_3d(species,candidate,flux,register_geometry, &
              interface_registers(level)%fields,level,accepted,boundary_faces,frame)
          else
            call synchronize_reactive_eb_tree_flux_3d(species,candidate,flux,accepted,level,level,boundary_faces,frame)
          endif
        endif
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) return
      else
      if(method/='none') before_sync=candidate
      if(old_only) then
        call move_alloc(physical_flux,sync_flux)
      else
        call move_alloc(flux,sync_flux)
      endif
      if(replayed_closure) then
        accepted=.true.
      else if(present(internal_energy).and.allocated(correction_work).and.method=='none') then
        if(observe_sync_enabled.and..not.allocated(diagnostic_raw)) allocate(diagnostic_raw,diagnostic_applied)
        if(observe_sync_enabled) call observe_mol_sync(0,level,time,interval, &
          candidate,held_candidate,sync_flux,scalar_flux)
        call close_held_prepared_eb_tree_3d(species,candidate,held_candidate,sync_flux,scalar_flux,held_policy, &
          correction_work,level,accepted,boundary_faces,frame, &
          extended_geometry=extended_geometry,internal_energy_change=clock_end_reset, &
          reset_current_level=godunov_mol.and.chemistry, &
          diagnostic_raw_increment=diagnostic_raw,diagnostic_applied_increment=diagnostic_applied)
        if(accepted.and.observe_sync_enabled) call observe_mol_sync(1,level,time,interval, &
          candidate,held_candidate,sync_flux,scalar_flux,diagnostic_raw,diagnostic_applied)
        paired_closed=accepted
      else
        call synchronize_reactive_eb_tree_flux_3d(species,candidate,sync_flux,accepted,level, &
          min(level+1,levels-1),boundary_faces,frame)
      endif
      if(.not.accepted.and.godunov_automatic_repair) go to 880
      if(.not.accepted.and.held_replay_enabled) then
        call repair_held_parent(history_first,sync_flux,scalar_flux,accepted)
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) then
          if(.not.held_replay_provider_refused) deferred_record(level)=history_first
          return
        endif
        physical_flux=sync_flux;replayed_closure=.true.
      endif
      if(.not.accepted.and.godunov_automatic_repair) go to 880
      if(.not.accepted.and.level+1<levels) then
        ! Paired ordinary held replay has already handled that family above;
        ! preserve the existing gas-only/nonreacting donor route here.
        if(present(internal_energy).and.chemistry) return
        call repair_reactive_eb_tree_interface_budget_3d(species,candidate,sync_flux,accepted,level,boundary_faces,frame)
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) return
        if(method/='none') before_sync=candidate
        call synchronize_reactive_eb_tree_flux_3d(species,candidate,sync_flux,accepted,level,level+1,boundary_faces,frame)
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
          correction%fields(p)%values=0
          target=tree%topology%nodes(p)%level
          if(target<level.or.target>level+1) cycle
          correction%fields(p)%values=candidate%nodes(p)%state-before_sync%nodes(p)%state
        end do
        ! This operator is linear in its argument, so an integrated reflux
        ! increment can be redistributed directly without division by dt.
        if(method=='state_redist') then
          call redistribute_eb_tree_state_3d(candidate,correction,redistributed,accepted,redist_target,boundary_faces, &
      max_order=redist_order)
        else
          call redistribute_sdc_eb_tree_flux_rate_3d(candidate,correction,redistributed,accepted,boundary_faces)
        end if
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) return
        candidate=before_sync
        do p=1,size(tree%nodes)
          target=tree%topology%nodes(p)%level
          if(target<level) then
            pending(target)%fields(p)%values=pending(target)%fields(p)%values+redistributed%fields(p)%values
          else
            candidate%nodes(p)%state=candidate%nodes(p)%state+redistributed%fields(p)%values+ &
              pending(level)%fields(p)%values
          end if
          pending(level)%fields(p)%values=0
        end do
        call recover_transfer(accepted)
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) return
        call average_down_reactive_eb_tree_3d(species,candidate,accepted,frame)
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) return
      end if
      endif
      if(present(internal_energy).and..not.paired_closed.and..not.replayed_closure) then
        if(prepared_redist.and.level+1<levels) then
          call synchronize_reactive_eb_tree_scalar_register_flux_3d(species,candidate,held_candidate,scalar_flux, &
            register_geometry,held_interface_registers(level)%fields,level,accepted,boundary_faces)
        else
          call synchronize_reactive_eb_tree_scalar_flux_3d(species,candidate,held_candidate,scalar_flux,accepted, &
            level,min(level+1,levels-1),boundary_faces)
        endif
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted.and.held_replay_enabled) then
          call repair_held_parent(history_first,flux,scalar_flux,accepted)
          if(.not.accepted.and.godunov_automatic_repair) go to 880
          if(.not.accepted) then
            if(.not.held_replay_provider_refused) deferred_record(level)=history_first
            return
          endif
          physical_flux=flux;replayed_closure=.true.
        endif
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) return
      endif
      if(allocated(correction_work).and..not.paired_closed.and. &
        .not.(replayed_closure.and.allocated(funding_context_work))) then
        call reset_held_eb_tree_family_3d(species,candidate,held_candidate,reset_difference,held_policy,level, &
          reset_change,accepted,frame,leaf_only=.true., &
          internal_energy_change=scalar_reset_change)
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted.and.held_replay_enabled) then
          call repair_held_parent(history_first,flux,scalar_flux,accepted)
          if(.not.accepted.and.godunov_automatic_repair) go to 880
          if(.not.accepted) then
            if(.not.held_replay_provider_refused) deferred_record(level)=history_first
            return
          endif
          physical_flux=flux;replayed_closure=.true.
          if(allocated(funding_context_work)) then
            reset_change=0;scalar_reset_change=0
          else
          call reset_held_eb_tree_family_3d(species,candidate,held_candidate,reset_difference,held_policy,level, &
            reset_change,accepted,frame,leaf_only=.true., &
          internal_energy_change=scalar_reset_change)
          endif
        endif
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) return
        clock_end_reset=clock_end_reset+scalar_reset_change
        correction_work=correction_work+reset_change
      endif
      if(routing) then
        if(held_auxiliary) then
          call apply_routed_held_source(level,accepted, &
            scalar_change=scalar_reset_change)
          if(.not.accepted.and.godunov_automatic_repair) go to 880
          if(.not.accepted) return
          clock_end_reset=clock_end_reset+scalar_reset_change
        else
        do p=1,size(tree%nodes)
          candidate%nodes(p)%state=candidate%nodes(p)%state+auxiliary_pending(level)%fields(p)%values
          auxiliary_pending(level)%fields(p)%values=0
        end do
        call recover_transfer(accepted)
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) return
        call average_down_reactive_eb_tree_3d(species,candidate,accepted,frame)
        if(.not.accepted.and.godunov_automatic_repair) go to 880
        if(.not.accepted) return
        endif
      endif
      if(chemistry.and.level+1<levels) then
        call average_down_reactive_eb_tree_field_3d(candidate,ir,accepted,minimum_depth=level)
        if(.not.accepted) return
      endif
      if(present(heat_release).and.level+1<levels) then
        call average_down_reactive_eb_tree_field_3d(candidate,heat,accepted,minimum_depth=level)
        if(.not.accepted) return
      endif
      if(present(mol_mode)) then
        if(.not.(mol_mode.and..not.chemistry.and.iterations>1)) physical_flux=flux
      else
        physical_flux=flux
      endif
      if(held_capture_enabled) then
        if(replayed_closure.and.allocated(funding_context_work)) then
          ! Native replay has already applied the replacement root EOF.
          clock_end_reset=held_history(history_first)%held_clock_end_change
        else if(replayed_closure.and.prepared_redist) then
          ! The nested accepted prepared replay has already closed its root.
          ! Retain that replacement exit instead of the discarded old trial.
          call compute_reactive_eb_held_projection(0._dp,held_history(history_first)%held_clock_end_change,0._dp, &
            clock_end_reset,scalar_reset_change,accepted)
          if(.not.accepted) return
          clock_end_reset=scalar_reset_change
        endif
        held_history(history_first)%held_clock_end_change=clock_end_reset
        held_history(history_first)%held_clock_events_closed=.true.
      endif
      if(godunov_capture) then
        godunov_records(godunov_index)%held_clock_end_change=clock_end_reset
        godunov_records(godunov_index)%held_clock_events_closed=.true.
        call close_godunov_clock(godunov_records(godunov_index),physical_flux,scalar_flux,accepted)
        if(.not.accepted) return
      endif
      if(allocated(funding_context_work).and..not.replayed_closure) then
        call close_eb_held_funding_clock_faces_3d(tree%topology,funding_geometry,nv,funding_local,level, &
          physical_flux,scalar_flux,funding_context_work,funding_next,accepted)
        if(.not.accepted) return
        call move_eb_held_funding_context_3d(funding_next,funding_context_work)
      endif
      good=.true.
      return
880   continue
      call repair_godunov_parent(godunov_index,godunov_cursor,level,flux,physical_flux,scalar_flux,accepted)
      good=accepted
    end subroutine advance_level
    subroutine apply_routed_held_source(level,good,scalar_change)
      integer, intent(in) :: level
      logical, intent(out) :: good
      real(dp), optional, intent(inout) :: scalar_change
      real(dp) :: source_reset_change
      type(amr_tree_state_3d) :: gas
      integer :: p
      source_reset_change=0
      allocate(gas%fields(size(tree%nodes)))
      do p=1,size(tree%nodes)
        gas%fields(p)%values=candidate%nodes(p)%state
      enddo
      call apply_eb_held_auxiliary_increment_3d(gas,held_candidate,auxiliary_pending(level),good)
      if(.not.good) return
      do p=1,size(tree%nodes)
        candidate%nodes(p)%state=gas%fields(p)%values
        auxiliary_pending(level)%fields(p)%values=0
      enddo
      if(allocated(correction_work)) then
        call reset_held_eb_tree_family_3d(species,candidate,held_candidate,reset_difference,held_policy,-1, &
          reset_change,good,frame,leaf_only=.true., &
          internal_energy_change=source_reset_change)
        if(.not.good) return
        correction_work=correction_work+reset_change
      else
        call recover_transfer(good)
        if(.not.good) return
      endif
      call average_down_reactive_eb_tree_3d(species,candidate,good,frame,held_candidate)
      if(.not.good) return
      if(present(scalar_change)) scalar_change=source_reset_change
    end subroutine apply_routed_held_source

    subroutine recover_transfer(good)
      logical, intent(out) :: good
      real(dp) :: primitive(reactive_nprim(size(species))),temperature,speed
      integer :: p,x,y,z,d(3)
      logical :: accepted
      good=.false.
      do p=1,size(tree%nodes)
        d=tree%topology%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(candidate%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
          if(any(candidate%nodes(p)%state(6:,x,y,z)<0._dp)) return
          if(held_auxiliary) then
          call recover_reactive_eb_cell_3d(species,candidate%nodes(p)%geometry,[x,y,z],candidate%nodes(p)%state(:,x,y,z), &
            candidate%nodes(p)%temperature(x,y,z),primitive,temperature,speed,accepted,frame, &
            held_candidate%fields(p)%values(1,x,y,z))
          else
          call recover_reactive_eb_cell_3d(species,candidate%nodes(p)%geometry,[x,y,z],candidate%nodes(p)%state(:,x,y,z), &
            candidate%nodes(p)%temperature(x,y,z),primitive,temperature,speed,accepted,frame)
          endif
          if(.not.accepted) return
          candidate%nodes(p)%temperature(x,y,z)=temperature
        end do
        end do
        end do
      end do
      good=.true.
    end subroutine recover_transfer
  end subroutine advance_reactive_eb_tree_sdc_subcycle_3d

  subroutine advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,iterations,chemistry, &
      transport_flags,force,rtol,atol,tree,extended_geometry,boundary_faces, &
      boundary_change,wall_change,chemical_change,source_change,minimum_theta,ok,reaction_source, &
      difmag,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,wall_temperature, &
      target_depth,ancestor_start,ancestor_end,fractions_start,fractions_end,source_start,source_end, &
      integrated_flux,old_predictor,new_predictor,redistribution,redistribution_transfer,target_volume_fraction, &
      auxiliary_source,auxiliary_change, &
      evaluate_reaction,mol_mode,mol_lax_friedrichs,embedded_wall_temperature,embedded_noslip,frame, &
      state_redist_max_order,state_redist_registers,evaluate_spatial_diffusion,spatial_diffusion_halo_width, &
      constant_coefficients,heat_release,registered_flux,internal_energy,internal_energy_flux, &
      ancestor_internal_energy_start,ancestor_internal_energy_end,held_state_redist_registers, &
      energy_policy,energy_correction,prepared_spatial_weight,held_interval_history,replay_interval,provider_refused, &
      spatial_diffusion_is_source,auxiliary_has_internal_energy,internal_energy_pressure_work_change,spatial_flux_weight, &
      godunov_interval_history,godunov_replay_interval,godunov_replay_partial, &
      godunov_source_recovery,godunov_source_recovery_ready,source_funding_stock,source_funding_impulse, &
      internal_energy_projection_change,internal_energy_closure_change,freeze_native_external_sources,native_spatial_weight, &
      native_source_prefix,native_source_prefix_ready,native_source_replay, &
      separate_momentum,consistent_velocity,paired_pressure,connected_directions,observe_mol,diagnostic_time,observe_mol_reaction)
    procedure(eb_mol_reaction_observer), optional :: observe_mol_reaction
    procedure(eb_mol_stage_observer), optional :: observe_mol
    real(dp), optional, intent(in) :: diagnostic_time
    logical, optional, intent(in) :: separate_momentum,consistent_velocity,paired_pressure,connected_directions
    logical :: selected_momentum,selected_velocity,selected_pressure,selected_connected
    type(eb_native_mol_source_prefix_3d), optional, intent(in) :: native_source_replay
    type(eb_native_mol_source_prefix_3d), optional, intent(inout) :: native_source_prefix
    logical, optional, intent(out) :: native_source_prefix_ready
    logical :: native_prefix_agreed,native_source_capture,native_prefix_scope
    integer :: native_prefix_observed
    logical, optional, intent(in) :: freeze_native_external_sources
    real(dp), optional, intent(in) :: native_spatial_weight
    real(dp) :: native_mol_weight
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
    integer :: naux,auxiliary_alloc_status
    logical, optional, intent(in) :: spatial_diffusion_is_source
    logical :: external_diffusion_source
    type(eb_mol_held_interval_3d), optional, intent(inout) :: held_interval_history
    type(eb_mol_held_interval_3d), optional, intent(in) :: replay_interval
    logical, optional, intent(out) :: provider_refused
    type(eb_mol_held_interval_3d) :: held_history_work
    type(dual_energy_policy), optional, intent(in) :: energy_policy
    real(dp), optional, intent(in) :: prepared_spatial_weight
    real(dp), optional, intent(inout) :: energy_correction
    type(dual_energy_policy) :: held_policy
    type(amr_tree_state_3d) :: reset_delta,predictor_reset
    real(dp) :: initial_reset,predictor_reset_change,final_reset,post_reset,reset_change,numerical_change,accepted_predictor_reset
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
    real(dp), optional, intent(in) :: embedded_wall_temperature
    logical, optional, intent(in) :: mol_mode,mol_lax_friedrichs,embedded_noslip
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    procedure(eb_sdc_spatial_diffusion_source), optional :: evaluate_spatial_diffusion
    integer, optional, intent(in) :: spatial_diffusion_halo_width
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: dt,force(3),rtol,atol
    integer, intent(in) :: iterations
    logical, intent(in) :: chemistry,transport_flags(4)
    type(eb_reactive_tree_3d), intent(inout) :: tree
    type(eb_geometry_3d), intent(in) :: extended_geometry(:)
    character(len=*), intent(in) :: boundary_faces(6)
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemical_change(:),source_change(:),minimum_theta
    ! Selected-level MOL stages publish their actual AMReX-style registers.
    ! The interval caller owns their later accumulation/restriction/reflux.
    type(state_redist_register_3d), allocatable, optional, intent(inout) :: state_redist_registers(:)
    type(state_redist_register_3d), allocatable :: stage_registers0(:),stage_registers1(:),final_registers(:)
    logical :: prepared_redist,prepared_godunov
    integer :: godunov_padding
    type(state_redist_register_3d), allocatable :: hydro_registers(:),iteration_registers(:),prior_hydro_registers(:)
    type(eb_geometry_3d), allocatable :: register_geometry(:)
    integer, optional, intent(in) :: state_redist_max_order
    logical, intent(out) :: ok
    integer, optional, intent(in) :: target_depth
    type(eb_reactive_tree_3d), optional, intent(in) :: ancestor_start,ancestor_end
    real(dp), optional, intent(in) :: fractions_start(:),fractions_end(:)
    type(amr_tree_state_3d), optional, intent(in) :: source_start,source_end
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: integrated_flux(:)
    type(amr_tree_state_3d), optional, intent(inout) :: old_predictor,new_predictor
    type(amr_tree_state_3d), optional, intent(inout) :: redistribution_transfer
    type(amr_tree_state_3d) :: nonface0,nonface1
    type(amr_tree_state_3d) :: transfer0,transfer1,hydro_transfer,final_transfer
    type(amr_tree_state_3d) :: admitted_mean_rate,admitted_mean_transfer
    type(eb_integrated_flux_3d), allocatable :: admitted_mean_faces(:),admitted_held_faces(:)
    type(amr_tree_state_3d) :: admitted_held_rate
    real(dp), optional, intent(in) :: target_volume_fraction
    procedure(eb_sdc_auxiliary_source), optional :: auxiliary_source
    real(dp), optional, intent(inout) :: auxiliary_change(:)
    type(amr_tree_state_3d) :: a0,a1
    real(dp), allocatable :: auxiliary(:),auxiliary_inventory0(:),auxiliary_inventory1(:)
    integer :: redist_order
    real(dp) :: redist_target
    character(len=*), optional, intent(in) :: redistribution
    character(len=16) :: redistribution_method
    type(amr_tree_state_3d), allocatable :: redist_halos
    type(eb_geometry_3d), allocatable :: redist_geometry(:)
    type(amr_tree_state_3d) :: redistributed,sync_rate
    type(eb_reactive_tree_3d) :: before_sync
    type(amr_tree_state_3d) :: op,np
    real(dp), allocatable :: midpoint(:)
    integer :: selected_depth
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source
    real(dp), optional, intent(in) :: difmag,inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: wall_temperature(6)
    type(eb_reactive_tree_3d) :: old,candidate,ledger
    type(amr_tree_state_3d) :: e0,e1,d0,d1,rates,hydro,ir,u,t
    type(amr_tree_state_3d) :: held_old,held_candidate,held_rate0,held_rate1,auxiliary_packed0,auxiliary_packed1
    type(eb_integrated_flux_3d), allocatable :: f0(:),f1(:),hf(:),baseline_hf(:)
    type(eb_integrated_flux_3d), allocatable :: held_flux0(:),held_flux1(:),held_flux(:)
    type(state_redist_register_3d), allocatable :: held_registers0(:),held_registers1(:),held_registers(:)
    real(dp) :: boundary(size(boundary_change)),wall(size(wall_change)),chemical(size(chemical_change))
    real(dp) :: external(size(source_change)),v0(size(source_change)),v1(size(source_change)),theta,stage_theta
    real(dp), allocatable :: primitive(:),wall_rate(:,:,:,:)
    real(dp) :: recovered,speed
    integer :: nv,node,n(3),it,i,j,k,c,sync_depth,callback_width
    logical :: valid,diffusive,mol,noslip,wall_active
    logical :: mol_source_funded,source_funded
    logical :: unsafe_forward,correction_periodic(3)
    real(dp) :: forward_state(reactive_nvar(size(species)))
    type(eb_geometry_3d), allocatable :: mol_geometry(:),heat_geometry(:)
    real(dp), allocatable :: heat_stage,heat_frequency,wall_impulse(:)
    real(dp) :: heat_old,heat_new,impulse_old(3),impulse_new(3)
    real(dp) :: mol_wall_old(size(wall_change)),mol_wall_new(size(wall_change))
    type(eb_mol_held_prepared_stage_3d), allocatable :: held_stages(:)
    type(eb_integrated_flux_3d), allocatable :: held_low_faces(:),held_low_scalar_faces(:)
    procedure(forced_reactor_integrator), pointer :: held_prepared_integrator
    real(dp) :: held_spatial_weight
    logical :: held_prepared,prepared_reacted,prepared_provider_refused,auxiliary_present
    godunov_partial=.false.
    if(present(godunov_replay_partial)) godunov_partial=godunov_replay_partial
    godunov_stage_capture=present(godunov_interval_history).or.present(godunov_source_recovery)
    godunov_recovery_agreed=.false.;godunov_recovery_forbidden=.false.;godunov_source_attempt_incomplete=.false.
    godunov_packet_pending=.false.
    if(present(godunov_source_recovery_ready)) godunov_source_recovery_ready=.false.
    native_prefix_agreed=.false.;native_prefix_observed=0
    native_source_capture=present(held_interval_history).or.present(native_source_prefix)
    if(present(native_source_prefix_ready)) native_source_prefix_ready=.false.
    ok=.false.
    call run_common()
    if(.not.ok.and.godunov_recovery_agreed) call publish_source_prefix()
    if(.not.ok.and.native_prefix_agreed) call publish_native_source_prefix()
  contains
    subroutine run_common()
    if(present(native_source_prefix).neqv.present(native_source_prefix_ready)) return
    native_prefix_agreed=present(native_source_prefix)
    if(present(godunov_source_recovery).neqv.present(godunov_source_recovery_ready)) return
    if(godunov_partial.and..not.present(godunov_replay_interval)) return
    godunov_recovery_agreed=present(godunov_source_recovery)
    held_auxiliary=.false.
    if(present(auxiliary_has_internal_energy)) held_auxiliary=auxiliary_has_internal_energy
    external_diffusion_source=.false.
    if(present(spatial_diffusion_is_source)) external_diffusion_source=spatial_diffusion_is_source
    if(external_diffusion_source.and..not.present(evaluate_spatial_diffusion)) then
      ok=.false.;return
    endif
    ok=.false.;nv=reactive_nvar(size(species))
    naux=eb_auxiliary_nvar(nv,auxiliary_has_internal_energy);theta=1
    allocate(auxiliary(naux) ,auxiliary_inventory0(naux),auxiliary_inventory1(naux),stat=auxiliary_alloc_status)
    if(auxiliary_alloc_status/=0) return
    if(held_auxiliary.and.(.not.present(internal_energy).or..not.present(auxiliary_source))) return
    if(present(provider_refused)) provider_refused=.false.
    held_spatial_weight=1;prepared_provider_refused=.false.
    native_external_frozen=.false.
    if(present(freeze_native_external_sources)) native_external_frozen=freeze_native_external_sources
    if(native_external_frozen.and..not.present(replay_interval)) return
    native_mol_weight=1
    if(present(native_spatial_weight)) then
      if(.not.native_external_frozen.and..not.present(native_source_replay)) return
      if(.not.ieee_is_finite(native_spatial_weight)) return
      if(native_spatial_weight<0.or.native_spatial_weight>1) return
      native_mol_weight=native_spatial_weight
    endif
    auxiliary_present=present(auxiliary_source)
    nullify(held_prepared_integrator)
    if(present(evaluate_reaction)) held_prepared_integrator=>classified_held_prepared_reactor
    initial_reset=0;predictor_reset_change=0;final_reset=0;post_reset=0;accepted_predictor_reset=0
    initial_scalar_reset=0
    mol_source_funded=.false.;source_funded=.false.
    held_policy=dual_energy_policy()
    if(present(energy_policy)) then
      if(.not.present(energy_correction)) return
      held_policy=energy_policy
    endif
    if(present(energy_correction)) then
      if(.not.present(internal_energy)) return
      if(.not.ieee_is_finite(energy_correction)) return
    endif
    if(present(constant_coefficients)) then
      if(.not.valid_constant_transport_coefficients(constant_coefficients)) return
    endif
    callback_width=6
    if(present(spatial_diffusion_halo_width)) callback_width=spatial_diffusion_halo_width
    if(callback_width<6.or.callback_width>huge(1)/4) return
    if(callback_width>6.and..not.present(evaluate_spatial_diffusion)) return
    prepared_redist=present(state_redist_registers)
    redist_order=0
    if(present(state_redist_max_order)) redist_order=state_redist_max_order
    if(redist_order/=0.and.redist_order/=2.and.redist_order/=4) return
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    mol=.false.
    if(present(mol_mode)) mol=mol_mode
    selected_momentum=.false.;selected_velocity=.false.;selected_pressure=.false.;selected_connected=.false.
    if(present(separate_momentum)) selected_momentum=separate_momentum
    if(present(consistent_velocity)) selected_velocity=consistent_velocity
    if(present(paired_pressure)) selected_pressure=paired_pressure
    if(present(connected_directions)) selected_connected=connected_directions
    if(selected_velocity.and..not.selected_momentum) return
    if(selected_connected.and..not.selected_momentum) return
    if(selected_pressure.and..not.selected_velocity) return
    if(selected_pressure.and.present(frame)) return
    if(present(mol_lax_friedrichs)) then
      if(mol_lax_friedrichs.and..not.mol) return
    end if
    if(present(frame)) then
      if(.not.mol) return
    end if
    mol_wall_old=0;mol_wall_new=0
    heat_old=0;heat_new=0;impulse_old=0;impulse_new=0
    noslip=.false.
    if(present(embedded_noslip)) noslip=embedded_noslip
    wall_active=noslip.or.present(embedded_wall_temperature)
    if(noslip.and..not.transport_flags(1)) return
    if(present(embedded_wall_temperature)) then
      if(.not.ieee_is_finite(embedded_wall_temperature).or.embedded_wall_temperature<=0) return
      if(.not.transport_flags(2)) return
    end if
    if(present(evaluate_reaction).and..not.chemistry) return
    if(.not.valid_tree(tree,nv)) return
    if(any([size(boundary_change),size(wall_change),size(chemical_change),size(source_change)]/=nv)) return
    if(iterations<1.or..not.all(ieee_is_finite([dt,force,rtol,atol]))) return
    if(min(dt,rtol,atol)<=0.or..not.valid_reactive_boundaries_3d(boundary_faces)) return
    if(transport_flags(4).and..not.transport_flags(3)) return
    if(present(auxiliary_source).neqv.present(auxiliary_change)) return
    if(present(auxiliary_change)) then
      if(size(auxiliary_change)/=naux) return
    end if
    selected_depth=-1
    if(present(target_depth)) then
      selected_depth=target_depth
      if(selected_depth<0.or.selected_depth>=tree%topology%level_count()) return
      if(.not.present(integrated_flux)) return
    end if
    if(present(source_funding_stock).neqv.present(source_funding_impulse)) return
    native_held_funding=selected_depth<0.or.present(source_funding_impulse)
    if(present(source_funding_stock)) then
      if(.not.present(internal_energy)) return
      if(.not.source_funding_stock%is_valid(tree%topology)) return
      if(size(source_funding_stock%fields(1)%values,1)/=nv+1) return
      do node=1,size(tree%nodes)
        if(any(.not.ieee_is_finite(source_funding_stock%fields(node)%values))) return
        if(selected_depth<0.or.tree%topology%nodes(node)%level==selected_depth) then
          if(any(source_funding_stock%fields(node)%values/=0)) return
        endif
        ! Pending loans belong to exposed fluid inventory. Virtual covered
        ! cells and solid EB cells cannot own a debt that admission later
        ! clears without settlement.
        n=tree%topology%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(covered(tree%topology,node,[i,j,k]).or.tree%nodes(node)%geometry%volume_fraction(i,j,k)<=0) then
            if(any(source_funding_stock%fields(node)%values(:,i,j,k)/=0)) return
          endif
        enddo
        enddo
        enddo
      enddo
    endif
    redistribution_method='none'
    if(present(redistribution)) redistribution_method=redistribution
    if(selected_momentum) then
      if(.not.mol.or..not.present(internal_energy).or.redistribution_method/='none'.or. &
        present(state_redist_registers)) return
    endif
    if(redistribution_method/='none'.and.redistribution_method/='flux_redist'.and. &
      redistribution_method/='state_redist') return
    redist_target=.5_dp
    if(present(target_volume_fraction)) redist_target=target_volume_fraction
    if(.not.ieee_is_finite(redist_target)) return
    if(redist_target<=0.or.redist_target>1) return
    if(selected_depth<0.and.mol.and.redistribution_method=='state_redist') prepared_redist=.true.
    if(selected_depth<0.and..not.mol.and.redistribution_method=='state_redist') then
      if(size(extended_geometry)==size(tree%nodes)) then
        prepared_godunov=.true.
        do node=1,size(tree%nodes)
          n=tree%topology%nodes(node)%shape
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
      godunov_padding=(extended_geometry(1)%nx-tree%topology%nodes(1)%shape(1))/2
    endif
    if(prepared_redist) then
      if(redistribution_method/='state_redist') return
      if(present(redistribution_transfer)) return
    else if(selected_depth>=0.and.redistribution_method/='none') then
      if(.not.present(redistribution_transfer)) return
    else
      if(present(redistribution_transfer)) return
    end if
    if(present(ancestor_start).or.present(ancestor_end).or.present(fractions_start).or.present(fractions_end)) then
      if(selected_depth<0) return
      if(.not.(present(ancestor_start).and.present(ancestor_end).and. &
        present(fractions_start).and.present(fractions_end))) return
      if(size(fractions_start)/=selected_depth.or.size(fractions_end)/=selected_depth) return
      if(.not.all(ieee_is_finite(fractions_start)).or..not.all(ieee_is_finite(fractions_end))) return
      if(any(fractions_start<0).or.any(fractions_end>1).or.any(fractions_end<fractions_start)) return
    end if
    if(present(source_start).neqv.present(source_end)) return
    if(present(source_start)) then
      if(.not.present(fractions_start).or..not.present(fractions_end)) return
      midpoint=.5_dp*(fractions_start+fractions_end)
    end if
    if(present(internal_energy_pressure_work_change).and..not.present(internal_energy)) return
    if(present(internal_energy_projection_change).and..not.present(internal_energy)) return
    if(present(internal_energy_closure_change).and..not.present(internal_energy)) return
    if(present(internal_energy).neqv.present(internal_energy_flux)) return
    if(present(held_state_redist_registers)) then
      if(.not.present(internal_energy).or.selected_depth<0.or..not.prepared_redist) return
    endif
    if(present(ancestor_internal_energy_start).neqv.present(ancestor_internal_energy_end)) return
    if(present(internal_energy)) then
      if(present(fractions_start).neqv.present(ancestor_internal_energy_start)) return
    endif
    if(present(ancestor_internal_energy_start)) then
      if(.not.present(internal_energy).or.selected_depth<0.or..not.present(fractions_start)) return
      if(.not.ancestor_internal_energy_start%is_valid(tree%topology).or. &
        .not.ancestor_internal_energy_end%is_valid(tree%topology)) return
      if(size(ancestor_internal_energy_start%fields(1)%values,1)/=1.or. &
        size(ancestor_internal_energy_end%fields(1)%values,1)/=1) return
    endif
    if(present(internal_energy)) then
      if((chemistry.and..not.present(energy_correction)).or. &
        (.not.mol.and.redistribution_method=='state_redist'.and..not.prepared_godunov).or. &
        (redistribution_method/='none'.and.redistribution_method/='state_redist').or. &
        (selected_depth>=0.and.redistribution_method=='state_redist'.and. &
        .not.present(held_state_redist_registers)).or. &
        (present(auxiliary_source).and..not.chemistry.and..not.held_auxiliary)) return
      if(.not.internal_energy%is_valid(tree%topology)) return
      if(size(internal_energy%fields(1)%values,1)/=1) return
      do node=1,size(tree%nodes)
        if(.not.all(ieee_is_finite(internal_energy%fields(node)%values))) return
      enddo
      held_old=internal_energy
      held_candidate=held_old
    endif
    if(present(held_interval_history)) then
      if(.not.mol.or..not.chemistry.or..not.present(internal_energy).or. &
        .not.present(energy_correction).or.selected_depth<0) return
      if(redistribution_method/='none'.and..not.prepared_redist) return
    endif
    if(godunov_stage_capture) then
      if(mol.or..not.present(internal_energy)) return
      if(redistribution_method/='none'.and..not.prepared_godunov) return
      if(present(replay_interval).or.present(held_interval_history)) return
      if(present(ancestor_start).and. &
        (.not.present(ancestor_internal_energy_start).or..not.present(ancestor_internal_energy_end))) return
    endif
    if(present(spatial_flux_weight)) then
      if(mol.or.present(prepared_spatial_weight).or.present(replay_interval)) return
      if(redistribution_method/='none'.and..not.prepared_godunov) return
      if(.not.ieee_is_finite(spatial_flux_weight)) return
      if(spatial_flux_weight<0.or.spatial_flux_weight>1) return
    endif
    native_prefix_scope=native_prefix_agreed.or.present(native_source_replay)
    if(native_prefix_scope) then
      ! First scope: selected ordinary clocks with built-in spatial physics.
      if(.not.mol.or..not.chemistry.or.selected_depth<0.or.redistribution_method/='none') return
      if(.not.present(internal_energy).or..not.present(energy_correction)) return
      if(present(evaluate_spatial_diffusion).or.present(godunov_source_recovery)) return
    endif
    held_prepared=mol.and.prepared_redist.and.chemistry.and.present(internal_energy)
    if(native_external_frozen.and.held_prepared) return
    if(present(prepared_spatial_weight)) then
      if(.not.held_prepared) return
      if(.not.ieee_is_finite(prepared_spatial_weight)) return
      if(prepared_spatial_weight<0.or.prepared_spatial_weight>1) return
      held_spatial_weight=prepared_spatial_weight
    endif
    if(present(replay_interval)) then
      if(.not.mol.or..not.chemistry.or..not.present(internal_energy).or. &
        .not.present(energy_correction).or.selected_depth<0) return
      if(replay_interval%spatial_diffusion_is_source.neqv.external_diffusion_source) return
      if(replay_interval%auxiliary_has_internal_energy.neqv.held_auxiliary) return
      if(held_prepared) then
        if(.not.present(prepared_spatial_weight)) return
        call validate_held_prepared_replay(valid)
      else
        if(redistribution_method/='none'.or.present(prepared_spatial_weight)) return
        call validate_held_native_replay(valid)
      endif
      if(.not.valid) return
    endif
    if(held_prepared) allocate(held_stages(0:iterations))
    if(present(godunov_replay_interval)) then
      if(mol.or..not.present(internal_energy).or.present(held_interval_history).or.present(replay_interval)) return
      if(redistribution_method/='none'.and..not.prepared_godunov) return
      call validate_godunov_replay(valid)
      if(.not.valid) return
    endif
    old=tree
    if(present(energy_correction)) then
      call reset_held_eb_tree_family_3d(species,old,held_old,reset_delta,held_policy,selected_depth, &
        initial_reset,valid,frame, &
        internal_energy_change=initial_scalar_reset,reset_level=merge(-2,-1,mol.and.chemistry.and.redistribution_method=='none'))
      if(.not.valid) return
    endif
    if(selected_depth<0.and..not.(mol.and.chemistry.and.redistribution_method=='none'.and.present(energy_correction))) then
      if(present(internal_energy)) then
        call average_down_reactive_eb_tree_3d(species,old,valid,frame,held_old)
      else
        call average_down_reactive_eb_tree_3d(species,old,valid,frame)
      endif
      if(.not.valid) return
    end if
    if(present(internal_energy)) held_candidate=held_old
    if(prepared_redist) then
      if(size(extended_geometry)/=size(tree%nodes)) return
      allocate(register_geometry(size(tree%nodes)))
      do node=1,size(tree%nodes)
        n=tree%topology%nodes(node)%shape
        if(prepared_godunov) then
          if(any([extended_geometry(node)%nx,extended_geometry(node)%ny,extended_geometry(node)%nz]/=n+2*godunov_padding)) return
          call extract_eb_geometry_3d(extended_geometry(node),[1,1,1]*(godunov_padding-5), &
            n+godunov_padding+6,register_geometry(node),valid)
          if(.not.valid) return
        else
          register_geometry(node)=extended_geometry(node)
        endif
      enddo
    endif
    if(mol.or.wall_active) then
      if(size(extended_geometry)/=size(tree%nodes)) return
      if(mol) allocate(mol_geometry(size(tree%nodes)))
      if(wall_active) then
        allocate(heat_geometry(size(tree%nodes)),heat_stage,heat_frequency)
        if(noslip) allocate(wall_impulse(3))
      end if
      do node=1,size(tree%nodes)
        n=tree%topology%nodes(node)%shape
        if(any([extended_geometry(node)%nx,extended_geometry(node)%ny,extended_geometry(node)%nz]/= &
          n+2*merge(godunov_padding,6,prepared_godunov))) return
        if(mol) then
          call extract_eb_geometry_3d(extended_geometry(node),[4,4,4],n+9,mol_geometry(node),valid)
          if(.not.valid) return
        end if
        if(wall_active) then
          if(prepared_godunov) then
            call extract_eb_geometry_3d(register_geometry(node),[5,5,5],n+8,heat_geometry(node),valid)
          else
            call extract_eb_geometry_3d(extended_geometry(node),[5,5,5],n+8,heat_geometry(node),valid)
          endif
          if(.not.valid) return
        end if
      end do
    end if
    allocate(primitive(reactive_nprim(size(species))),ir%fields(size(tree%nodes)))
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      allocate(ir%fields(node)%values(nv-4,n(1),n(2),n(3)),source=0._dp)
    end do
    if(present(heat_release)) then
      if(.not.heat_release%is_valid(tree%topology)) return
      if(size(heat_release%fields(1)%values,1)/=1) return
      allocate(heat);heat=heat_release
    endif
    if(present(heat_release).and..not.chemistry) then
      do node=1,size(tree%topology%nodes)
        if(selected_depth>=0.and.tree%topology%nodes(node)%level/=selected_depth) cycle
        heat%fields(node)%values=0
      enddo
    endif
    if(present(reaction_source)) then
      if(.not.allocated(reaction_source%fields)) return
      if(size(reaction_source%fields)/=size(tree%nodes)) return
      do node=1,size(tree%nodes)
        if(.not.allocated(reaction_source%fields(node)%values)) return
        if(any(shape(reaction_source%fields(node)%values)/=shape(ir%fields(node)%values))) return
      end do
      ir=reaction_source
    end if
    do node=1,size(tree%nodes)
      if(.not.active(node)) cycle
      ! The pinned disabled-chemistry path resets I_R before predicting.
      if(.not.chemistry) ir%fields(node)%values=0
      do c=1,nv-4
        where(old%nodes(node)%geometry%volume_fraction<=0) ir%fields(node)%values(c,:,:,:)=0
      end do
      if(.not.all(ieee_is_finite(ir%fields(node)%values))) return
    end do
    if(present(replay_interval)) then
      do node=1,size(tree%nodes)
        if(active(node)) ir%fields(node)=replay_interval%initial_reaction_source%fields(node)
      enddo
    endif
    if(present(native_source_replay)) then
      call validate_native_source_replay(valid)
      if(.not.valid) return
    endif
    if(native_source_capture) then
      held_history_work%prepared_redistribution=held_prepared
      held_history_work%accepted_coarse_stages=held_prepared
      held_history_work%spatial_diffusion_is_source=external_diffusion_source
      held_history_work%auxiliary_has_internal_energy=held_auxiliary
      held_history_work%replayed=present(replay_interval).or.present(native_source_replay)
      held_history_work%has_auxiliary_source=present(auxiliary_source)
      if(.not.held_prepared) then
        held_history_work%separate_momentum=selected_momentum
        held_history_work%consistent_velocity=selected_velocity
        held_history_work%paired_pressure=selected_pressure
        held_history_work%connected_directions=selected_connected
      held_history_work%native_iterations=iterations;held_history_work%native_force=force
        held_history_work%native_spatial_weight=native_mol_weight
        allocate(held_history_work%external_stages(0:iterations),stat=history_status)
        if(history_status/=0) return
        held_history_work%native_providers=[present(evaluate_spatial_diffusion),.false.]
        held_history_work%native_provider_halo_width=callback_width
        if(external_diffusion_source.and.any(held_history_work%native_providers)) then
          allocate(held_history_work%provider_stages(2,0:iterations),stat=history_status)
          if(history_status/=0) return
        endif
        if(present(auxiliary_source)) then
          allocate(held_history_work%auxiliary_stages(0:iterations),stat=history_status)
          if(history_status/=0) return
        endif
      endif
      held_history_work%initial=old;held_history_work%initial_internal_energy=held_old
      held_history_work%initial_reaction_source=ir;held_history_work%policy=held_policy
      held_history_work%clock_initial=old;held_history_work%clock_internal_energy=held_old
      held_history_work%clock_reaction_source=ir
      held_history_work%clock_energy_correction=energy_correction+initial_reset
      if(present(heat_release)) held_history_work%clock_heat=heat
      held_history_work%has_frame=present(frame)
      if(present(frame)) held_history_work%frame=frame
      if(present(ancestor_start)) then
        held_history_work%ancestor_start=ancestor_start;held_history_work%ancestor_end=ancestor_end
        held_history_work%fractions_start=fractions_start;held_history_work%fractions_end=fractions_end
        if(present(ancestor_internal_energy_start)) then
          held_history_work%ancestor_internal_energy_start=ancestor_internal_energy_start
          held_history_work%ancestor_internal_energy_end=ancestor_internal_energy_end
        endif
        if(present(source_start)) then
          held_history_work%ancestor_source_start=source_start;held_history_work%ancestor_source_end=source_end
        endif
      endif
    endif
    if(godunov_stage_capture) then
      call begin_godunov_history(valid)
      if(.not.valid) return
    endif
    candidate=old;diffusive=any(transport_flags(:3)).or.present(evaluate_spatial_diffusion)
    call sources(old,e0,d0,a0,0,f0,transfer0,valid,fractions_start)
    if(.not.valid) then
      godunov_source_attempt_incomplete=.true.;return
    endif
    if(native_prefix_agreed) native_prefix_observed=max(native_prefix_observed,1)
    if(godunov_stage_capture) then
      call capture_godunov_source(0,old,e0,d0,a0,f0,valid)
      if(.not.valid) then
        godunov_source_attempt_incomplete=.true.;return
      endif
    endif
    e1=e0;d1=d0;mol_wall_new=mol_wall_old;heat_new=heat_old;impulse_new=impulse_old
    if(prepared_redist) stage_registers1=stage_registers0
    if(present(auxiliary_source)) a1=a0
    call predictor_field(e0,d0,a0,op)
    allocate(u%fields(size(tree%nodes)),t%fields(size(tree%nodes)))
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      u%fields(node)%values=old%nodes(node)%state
      allocate(t%fields(node)%values(1,n(1),n(2),n(3)))
      t%fields(node)%values(1,:,:,:)=old%nodes(node)%temperature
    end do
    if(mol.and.present(internal_energy)) then
      held_hydro=held_old;held_hydro_flux=held_flux0
      do node=1,size(tree%nodes)
        held_hydro%fields(node)%values=0
        if(allocated(held_hydro_flux(node)%x)) held_hydro_flux(node)%x=0
        if(allocated(held_hydro_flux(node)%y)) held_hydro_flux(node)%y=0
        if(allocated(held_hydro_flux(node)%z)) held_hydro_flux(node)%z=0
      enddo
    endif
    if(mol) then
      hydro=d0;hf=f0
      do node=1,size(tree%nodes)
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
        do node=1,size(tree%nodes)
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
      ! Record reset forcing only for the final accepted SDC iterate.
      accepted_predictor_reset=0
      if(allocated(admitted_mean_rate%fields)) deallocate(admitted_mean_rate%fields)
      if(allocated(admitted_mean_transfer%fields)) deallocate(admitted_mean_transfer%fields)
      if(allocated(admitted_mean_faces)) deallocate(admitted_mean_faces)
      if(allocated(admitted_held_faces)) deallocate(admitted_held_faces)
      if(allocated(admitted_held_rate%fields)) deallocate(admitted_held_rate%fields)
      ! MOL builds an Euler predictor only for the first iteration. Later
      ! reacting iterations use the preceding reacted endpoint directly.
      if(.not.mol.or.it==1) then
      call combined_rates(.true.)
      if(.not.mol) then
      if(prepared_godunov) then
        if(it>1) prior_hydro_registers=hydro_registers
        if(it>1.and.present(internal_energy)) prior_held_hydro_registers=held_hydro_registers
        if(present(internal_energy)) then
      call compute_godunov_state_redist_eb_tree_hydro_3d(species,old,rates,extended_geometry,dt,hydro,hf,valid, &
          difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          target_depth,ancestor_start,ancestor_end,fractions_start,source_start,source_end,midpoint, &
          registers=iteration_registers,max_order=redist_order,target_volume_fraction=redist_target, &
          internal_energy_state=held_old,internal_energy_rate=held_hydro, &
          internal_energy_flux=held_hydro_flux,internal_energy_pressure_work=held_hydro_work, &
          ancestor_internal_energy_start=ancestor_internal_energy_start, &
          ancestor_internal_energy_end=ancestor_internal_energy_end, &
          internal_energy_registers=iteration_held_hydro_registers, &
        spatial_flux_weight=spatial_flux_weight)
      else
      call compute_godunov_state_redist_eb_tree_hydro_3d(species,old,rates,extended_geometry,dt,hydro,hf,valid, &
          difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          target_depth,ancestor_start,ancestor_end,fractions_start,source_start,source_end,midpoint, &
          registers=iteration_registers,max_order=redist_order,target_volume_fraction=redist_target, &
        spatial_flux_weight=spatial_flux_weight)
      endif
        if(.not.valid) return
        if(it>1) then
          do node=1,size(tree%nodes)
            if(.not.active(node)) cycle
            n=tree%topology%nodes(node)%shape
            call accumulate_state_redist_coarse_register_3d(iteration_registers(node),prior_hydro_registers(node),nv,n+12,valid)
            if(.not.valid) return
            if(present(internal_energy)) then
              call accumulate_state_redist_coarse_register_3d(iteration_held_hydro_registers(node), &
                prior_held_hydro_registers(node),1,n+12,valid)
              if(.not.valid) return
            endif
          enddo
        endif
        call move_alloc(iteration_registers,hydro_registers)
        if(present(internal_energy)) call move_alloc(iteration_held_hydro_registers,held_hydro_registers)
      else
      if(present(internal_energy)) then
      call compute_sdc_eb_tree_hydro_3d(species,old,rates,extended_geometry,dt,hydro,hf,valid, &
        difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        target_depth,ancestor_start,ancestor_end,fractions_start,source_start,source_end,midpoint, &
          internal_energy_state=held_old,internal_energy_rate=held_hydro, &
          internal_energy_flux=held_hydro_flux,internal_energy_pressure_work=held_hydro_work, &
          ancestor_internal_energy_start=ancestor_internal_energy_start, &
          ancestor_internal_energy_end=ancestor_internal_energy_end, &
        spatial_flux_weight=spatial_flux_weight)
      else
      call compute_sdc_eb_tree_hydro_3d(species,old,rates,extended_geometry,dt,hydro,hf,valid, &
        difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        target_depth,ancestor_start,ancestor_end,fractions_start,source_start,source_end,midpoint, &
        spatial_flux_weight=spatial_flux_weight)
      endif
      endif
      if(.not.valid) return
      unsafe_forward=.false.
      do node=1,size(tree%nodes)
        if(.not.active(node)) cycle
        n=tree%topology%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(old%nodes(node)%geometry%volume_fraction(i,j,k)<=0) cycle
          forward_state=old%nodes(node)%state(:,i,j,k)+ &
            dt*(rates%fields(node)%values(:,i,j,k)+hydro%fields(node)%values(:,i,j,k))
          if(any(.not.ieee_is_finite(forward_state))) return
          unsafe_forward=unsafe_forward.or.forward_state(1)<=0.or.any(forward_state(6:)<0)
        enddo
        enddo
        enddo
      enddo
      if(unsafe_forward.and.present(internal_energy).and.native_held_funding.and..not.prepared_redist) then
        call fund_held_native_stage(valid)
        if(valid) unsafe_forward=.false.
      endif
      if(unsafe_forward) then
        if(prepared_godunov) then
        if(present(internal_energy)) then
      call compute_godunov_state_redist_eb_tree_hydro_3d(species,old,rates,extended_geometry,dt,hydro,hf,valid, &
          difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          target_depth,ancestor_start,ancestor_end,fractions_start,source_start,source_end,midpoint, &
          registers=iteration_registers,max_order=redist_order,target_volume_fraction=redist_target,correct_forward_flux=.true., &
          internal_energy_state=held_old,internal_energy_rate=held_hydro, &
          internal_energy_flux=held_hydro_flux,internal_energy_pressure_work=held_hydro_work, &
          ancestor_internal_energy_start=ancestor_internal_energy_start, &
          ancestor_internal_energy_end=ancestor_internal_energy_end, &
          internal_energy_registers=iteration_held_hydro_registers, &
        spatial_flux_weight=spatial_flux_weight)
      else
      call compute_godunov_state_redist_eb_tree_hydro_3d(species,old,rates,extended_geometry,dt,hydro,hf,valid, &
          difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          target_depth,ancestor_start,ancestor_end,fractions_start,source_start,source_end,midpoint, &
          registers=iteration_registers,max_order=redist_order,target_volume_fraction=redist_target,correct_forward_flux=.true., &
        spatial_flux_weight=spatial_flux_weight)
      endif
        if(.not.valid) return
        if(it>1) then
          do node=1,size(tree%nodes)
            if(.not.active(node)) cycle
            n=tree%topology%nodes(node)%shape
            call accumulate_state_redist_coarse_register_3d(iteration_registers(node),prior_hydro_registers(node),nv,n+12,valid)
            if(.not.valid) return
            if(present(internal_energy)) then
              call accumulate_state_redist_coarse_register_3d(iteration_held_hydro_registers(node), &
                prior_held_hydro_registers(node),1,n+12,valid)
              if(.not.valid) return
            endif
          enddo
        endif
        call move_alloc(iteration_registers,hydro_registers)
        if(present(internal_energy)) call move_alloc(iteration_held_hydro_registers,held_hydro_registers)
      valid=.true.
      do node=1,size(tree%nodes)
        if(.not.active(node)) cycle
        n=tree%topology%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(old%nodes(node)%geometry%volume_fraction(i,j,k)<=0) cycle
          forward_state=old%nodes(node)%state(:,i,j,k)+ &
            dt*(rates%fields(node)%values(:,i,j,k)+hydro%fields(node)%values(:,i,j,k))
          if(any(.not.ieee_is_finite(forward_state))) return
          valid=valid.and.forward_state(1)>0.and.all(forward_state(6:)>=0)
        enddo
        enddo
        enddo
      enddo
      if(.not.valid) return
        else
        if(present(internal_energy)) then
      call compute_sdc_eb_tree_hydro_3d(species,old,rates,extended_geometry,dt,hydro,hf,valid, &
          difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          target_depth,ancestor_start,ancestor_end,fractions_start,source_start,source_end,midpoint,baseline_hf, &
          internal_energy_state=held_old,internal_energy_rate=held_hydro, &
          internal_energy_flux=held_hydro_flux,internal_energy_pressure_work=held_hydro_work, &
          ancestor_internal_energy_start=ancestor_internal_energy_start, &
          ancestor_internal_energy_end=ancestor_internal_energy_end,internal_energy_baseline_flux=held_hydro_baseline, &
        spatial_flux_weight=spatial_flux_weight)
      else
      call compute_sdc_eb_tree_hydro_3d(species,old,rates,extended_geometry,dt,hydro,hf,valid, &
          difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          target_depth,ancestor_start,ancestor_end,fractions_start,source_start,source_end,midpoint,baseline_hf, &
        spatial_flux_weight=spatial_flux_weight)
      endif
        if(.not.valid) return
        correction_periodic=boundary_faces([1,3,5])=='periodic'.and.boundary_faces([2,4,6])=='periodic'
        if(present(internal_energy)) then
        call apply_eb_tree_flux_correction_3d(old,rates,baseline_hf,dt,hydro,hf,valid,correction_periodic,target_depth, &
          internal_energy_baseline=held_hydro_baseline,internal_energy_rate=held_hydro,internal_energy_flux=held_hydro_flux)
        else
        call apply_eb_tree_flux_correction_3d(old,rates,baseline_hf,dt,hydro,hf,valid,correction_periodic,target_depth)
        endif
        if(.not.valid) return
        endif
      endif
      if(redistribution_method/='none'.and..not.prepared_godunov) then
        if(redistribution_method=='state_redist'.and.redist_order>0) then
          call redist_boundary_samples(old,valid)
          if(.not.valid) return
        endif
        if(selected_depth>=0) then
          if(redistribution_method=='state_redist') then
            call redistribute_sdc_eb_tree_state_level_rate_3d(old,hydro,selected_depth,dt,redistributed, &
              hydro_transfer,valid,redist_target,boundary_faces, &
        max_order=redist_order,boundary_halos=redist_halos,boundary_geometry=redist_geometry)
          else
            call redistribute_sdc_eb_tree_level_rate_3d(old,hydro,selected_depth,dt,redistributed, &
              hydro_transfer,valid,boundary_faces)
          end if
        else
          if(redistribution_method=='state_redist') then
            call redistribute_sdc_eb_tree_state_rate_3d(old,hydro,dt,redistributed,valid,redist_target,boundary_faces, &
        max_order=redist_order,boundary_halos=redist_halos,boundary_geometry=redist_geometry)
          else
            call redistribute_sdc_eb_tree_flux_rate_3d(old,hydro,redistributed,valid,boundary_faces)
          end if
        end if
        if(.not.valid) return
        hydro=redistributed
      end if
      ! EB.cpp posts spatial density even when no redistribution is requested.
      do node=1,size(tree%nodes)
        if(.not.active(node)) cycle
        call post_eb_spatial_density_rate_3d(old%nodes(node)%geometry,hydro%fields(node)%values,valid)
        if(.not.valid) return
      enddo
      end if
      if(mol.and.present(internal_energy).and.native_held_funding.and..not.prepared_redist) then
        unsafe_forward=.false.
        do node=1,size(tree%nodes)
          if(.not.active(node)) cycle
          n=tree%topology%nodes(node)%shape
          do k=1,n(3)
          do j=1,n(2)
          do i=1,n(1)
            if(old%nodes(node)%geometry%volume_fraction(i,j,k)<=0) cycle
            forward_state=old%nodes(node)%state(:,i,j,k)+ &
              dt*(rates%fields(node)%values(:,i,j,k)+hydro%fields(node)%values(:,i,j,k))
            unsafe_forward=unsafe_forward.or.any(.not.ieee_is_finite(forward_state)).or. &
              forward_state(1)<=0.or.any(forward_state(6:)<0)
          enddo
          enddo
          enddo
        enddo
        if(unsafe_forward) then
          call fund_held_native_stage(valid)
          if(.not.valid) return
        endif
      endif
      if(.not.mol.and.present(internal_energy)) then
        godunov_predictor_forcing=rates;godunov_predictor_hydro=hydro
        godunov_predictor_held=held_hydro;godunov_predictor_work=held_hydro_work
        godunov_predictor_faces=hf;godunov_predictor_scalar_faces=held_hydro_flux
      endif
      if(godunov_stage_capture) then
        call capture_godunov_hydro(it,valid)
        if(.not.valid) return
      endif
      do node=1,size(tree%nodes)
        if(.not.active(node)) cycle
        candidate%nodes(node)%state=old%nodes(node)%state+dt*(rates%fields(node)%values+hydro%fields(node)%values)
      end do
      if(present(energy_correction)) then
        call held_update(held_old,held_rate0,held_rate0,.5_dp,held_candidate,valid)
        if(.not.valid) return
      endif
      call recover(valid,1)
      if(.not.valid) return
      end if
      if(present(internal_energy).and..not.present(energy_correction)) then
        call held_update(held_old,held_rate0,held_rate0,.5_dp,held_candidate,valid)
        if(.not.valid) return
      endif
      call sources(candidate,e1,d1,a1,it,f1,transfer1,valid,fractions_end)
      if(.not.valid) then
        godunov_source_attempt_incomplete=.true.;return
      endif
      if(native_prefix_agreed) native_prefix_observed=max(native_prefix_observed,it+1)
      if(godunov_stage_capture) then
        call capture_godunov_source(it,candidate,e1,d1,a1,f1,valid)
        if(.not.valid) then
          godunov_source_attempt_incomplete=.true.;return
        endif
      endif
      if(mol.and.it>1) then
        ! Pinned do_mol_advance rebuilds only the spatial new-time RHS on
        ! later reacting passes. molSrc_old retains old external forcing;
        ! the first predictor's new external forcing is not added again.
        do node=1,size(tree%nodes)
          if(.not.active(node)) cycle
          e1%fields(node)%values=0
        enddo
      endif
      call combined_rates(.false.)
      prepared_reacted=.false.
      if(held_prepared) then
        if(present(replay_interval)) then
          ! A clock-group trial owns this weight. Never repair an individual
          ! failed stage with a different family or retry an admitted refusal.
          call run_held_prepared_reaction(valid)
          prepared_reacted=valid
          theta=min(theta,held_spatial_weight)
        else
          call limit_held_prepared_pass(it,valid,prepared_reacted)
        endif
        if(.not.valid) return
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
          if(present(energy_correction).and.it==1.and..not.held_prepared) then
            do node=1,size(tree%nodes)
            if(.not.active(node)) cycle
              rates%fields(node)%values(5,:,:,:)=rates%fields(node)%values(5,:,:,:)+ &
                .5_dp*predictor_reset%fields(node)%values(5,:,:,:)/dt
            enddo
          endif
          call admit_held_native_corrector(valid)
          rates=unreset_rates
          if(.not.valid) return
          end block
        endif
        call limit_mol_pass(valid)
        if(.not.valid) return
        call combined_rates(.false.)
      endif
      if(present(energy_correction).and.it==1.and..not.held_prepared) then
        accepted_predictor_reset=.5_dp*predictor_reset_change
        do node=1,size(tree%nodes)
          if(.not.active(node)) cycle
          rates%fields(node)%values(5,:,:,:)=rates%fields(node)%values(5,:,:,:)+ &
            .5_dp*predictor_reset%fields(node)%values(5,:,:,:)/dt
        enddo
      endif
      if(present(internal_energy).and.native_held_funding.and..not.prepared_redist) then
        call admit_held_native_corrector(valid)
        if(.not.valid) return
      endif
      do node=1,size(tree%nodes)
        rates%fields(node)%values=rates%fields(node)%values+hydro%fields(node)%values
      end do
      if(chemistry) then
        if(.not.prepared_reacted) then
        if(present(internal_energy)) then
        if(native_held_funding.and..not.prepared_redist) then
          call prepare_covered_reactor_forcing(valid)
          if(.not.valid) return
        endif
        call react_sdc_eb_tree_3d(species,reactions,old,rates,dt,rtol,atol,u,t,ir,valid,target_depth, &
          evaluate_reaction=godunov_reactor_adapter,frame=frame,heat_release=heat, &
          internal_energy=held_candidate,old_internal_energy=held_old)
        else
        call react_sdc_eb_tree_3d(species,reactions,old,rates,dt,rtol,atol,u,t,ir,valid,target_depth, &
          evaluate_reaction=evaluate_reaction,frame=frame,heat_release=heat)
        endif
        if(.not.valid) return
        if(mol.and.present(observe_mol_reaction)) then
          block
            real(dp) :: observed_time
            observed_time=0
            if(present(diagnostic_time)) observed_time=diagnostic_time
            call observe_mol_reaction(it,selected_depth,observed_time,dt,old,u,rates,d0,d1,e0,e1,ir)
          end block
        endif
        do node=1,size(tree%nodes)
          if(.not.active(node)) cycle
          candidate%nodes(node)%state=u%fields(node)%values
          candidate%nodes(node)%temperature=t%fields(node)%values(1,:,:,:)
        end do
        if(present(energy_correction)) then
          call recover(valid,2)
          if(.not.valid) return
        endif
        endif
      else
        do node=1,size(tree%nodes)
          if(.not.active(node)) cycle
          candidate%nodes(node)%state=old%nodes(node)%state+dt*rates%fields(node)%values
        end do
        if(present(energy_correction)) then
          call held_update(held_old,held_rate0,held_rate1,.5_dp,held_candidate,valid)
          if(.not.valid) return
        endif
        call recover(valid,2)
        if(.not.valid) return
      end if
      if(present(internal_energy).and..not.present(energy_correction)) then
        call held_update(held_old,held_rate0,held_rate1,.5_dp,held_candidate,valid)
        if(.not.valid) return
      endif
    end do
    projection_inventory=0;closure_inventory=0
    if(present(internal_energy)) then
      call measure_held_projection(valid)
      if(.not.valid) return
    endif
    ledger=old
    do node=1,size(tree%nodes)
      ledger%nodes(node)%state=0
      if(.not.active(node)) cycle
      ledger%nodes(node)%state=candidate%nodes(node)%state-(old%nodes(node)%state+dt*rates%fields(node)%values)
    end do
    call integrate_reactive_eb_tree_3d(ledger,chemical,valid)
    if(.not.valid) return
    if(present(energy_correction)) chemical(5)=chemical(5)-final_reset
    call field_inventory(e0,v0,valid)
    if(.not.valid) return
    call field_inventory(e1,v1,valid)
    if(.not.valid) return
    external=.5_dp*dt*(v0+v1)
    if(present(auxiliary_source)) then
      if(held_auxiliary) then
        call field_inventory(auxiliary_packed0,auxiliary_inventory0,valid)
        if(.not.valid) return
        call field_inventory(auxiliary_packed1,auxiliary_inventory1,valid)
        if(.not.valid) return
        auxiliary=.5_dp*dt*(auxiliary_inventory0+auxiliary_inventory1)
      else
        call field_inventory(a0,v0,valid)
        if(.not.valid) return
        call field_inventory(a1,v1,valid)
        if(.not.valid) return
        auxiliary=.5_dp*dt*(v0+v1)
      endif
      if(.not.all(ieee_is_finite(auxiliary))) return
    end if
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      ledger%nodes(node)%state=0
      if(.not.active(node)) cycle
      allocate(wall_rate(nv,n(1),n(2),n(3)))
      call wall_source(old,node,wall_rate,valid,0)
      if(.not.valid) return
      ledger%nodes(node)%state=dt*wall_rate
      deallocate(wall_rate)
      if(diffusive.or.mol) then
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
    call integrate_reactive_eb_tree_3d(ledger,wall,valid)
    if(.not.valid) return
    if(mol) wall=.5_dp*(mol_wall_old+mol_wall_new)
    wall(5)=wall(5)+.5_dp*(heat_old+heat_new)
    wall(2:4)=wall(2:4)+.5_dp*(impulse_old+impulse_new)
    if(.not.all(ieee_is_finite(wall))) return
    if(present(registered_flux).or.godunov_stage_capture.or.(mol.and..not.chemistry.and.iterations>1)) then
      registration=hf
      if(mol.and..not.chemistry.and.iterations>1) then
        ! No later reacting pass exists to register the final new-time half.
        ! Assemble directly from f0, avoiding subtraction cancellation.
        do node=1,size(registration)
          if(.not.active(node)) cycle
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
      call ensure_selected_scalar_faces(valid)
      if(.not.valid) return
      held_flux=held_flux0
      do node=1,size(tree%nodes)
        if(.not.active(node).and..not.native_held_funding) cycle
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
        if(.not.all(ieee_is_finite(held_flux(node)%x)).or. &
          .not.all(ieee_is_finite(held_flux(node)%y)).or. &
          .not.all(ieee_is_finite(held_flux(node)%z))) return
      enddo
      if(prepared_redist) then
        if(prepared_godunov) then
          do node=1,size(held_registers0)
            if(allocated(held_registers0(node)%fine_delta)) &
              held_registers0(node)%fine_delta=.5_dp*held_registers0(node)%fine_delta
            if(allocated(held_registers1(node)%fine_delta)) &
              held_registers1(node)%fine_delta=.5_dp*held_registers1(node)%fine_delta
          enddo
        endif
        held_registers=held_registers0
        if(selected_depth>=0.and.mol.and..not.chemistry.and.iterations>1) then
          do node=1,size(held_registers1)
            if(allocated(held_registers1(node)%fine_delta)) held_registers1(node)%fine_delta=0
          enddo
        endif
        do node=1,size(tree%nodes)
          if(.not.active(node)) cycle
          n=tree%topology%nodes(node)%shape
          call accumulate_state_redist_register_3d(held_registers(node),held_registers1(node),1,n+12,valid)
          if(.not.valid) return
          if(prepared_godunov) then
            call accumulate_state_redist_register_3d(held_registers(node),held_hydro_registers(node),1,n+12,valid)
            if(.not.valid) return
          endif
        enddo
      endif
    endif
    if(prepared_redist) then
      final_registers=stage_registers0
      ! Selected-level registration follows the pinned zero new-time factor;
      ! a simultaneous hierarchy update retains both physical half steps.
      if(selected_depth>=0.and.mol.and..not.chemistry.and.iterations>1) then
        do node=1,size(stage_registers1)
          if(allocated(stage_registers1(node)%fine_delta)) stage_registers1(node)%fine_delta=0
        enddo
      endif
      if(prepared_godunov) then
        do node=1,size(tree%nodes)
          if(allocated(final_registers(node)%fine_delta)) final_registers(node)%fine_delta= &
            .5_dp*final_registers(node)%fine_delta
          if(allocated(stage_registers1(node)%fine_delta)) stage_registers1(node)%fine_delta= &
            .5_dp*stage_registers1(node)%fine_delta
        enddo
      endif
      do node=1,size(tree%nodes)
        if(.not.allocated(final_registers(node)%coarse_delta).and. &
          .not.allocated(final_registers(node)%fine_delta)) cycle
        n=tree%topology%nodes(node)%shape
        ! Match the accepted endpoint source in both coarse and fine ledgers.
        call accumulate_state_redist_register_3d(final_registers(node),stage_registers1(node),nv,n+12,valid)
        if(.not.valid) return
        if(prepared_godunov) then
          call accumulate_state_redist_register_3d(final_registers(node),hydro_registers(node),nv,n+12,valid)
          if(.not.valid) return
        endif
      enddo
    endif
    if(selected_depth<0.and.present(internal_energy).and.present(energy_correction)) then
      if(tree%topology%level_count()==1) then
        call close_held_prepared_eb_tree_3d(species,candidate,held_candidate,hf,held_flux,held_policy,post_reset, &
          0,valid,boundary_faces,frame, &
          extended_geometry=extended_geometry,internal_energy_change=closure_inventory(1), &
          reset_current_level=mol.and.chemistry.and.redistribution_method=='none')
        if(.not.valid) return
      endif
      do sync_depth=tree%topology%level_count()-2,0,-1
        if(prepared_redist) then
          call close_held_prepared_eb_tree_3d(species,candidate,held_candidate,hf,held_flux,held_policy,post_reset, &
            sync_depth,valid,boundary_faces,frame,register_geometry,final_registers,held_registers, &
          internal_energy_change=closure_inventory(1))
        else
          call close_held_prepared_eb_tree_3d(species,candidate,held_candidate,hf,held_flux,held_policy,post_reset, &
            sync_depth,valid,boundary_faces,frame, &
          extended_geometry=extended_geometry,internal_energy_change=closure_inventory(1), &
          reset_current_level=mol.and.chemistry.and.redistribution_method=='none')
        endif
        if(.not.valid) return
      enddo
      if(prepared_redist) then
      do node=1,size(tree%nodes)
        if(allocated(final_registers(node)%coarse_delta)) then
          if(any(final_registers(node)%coarse_delta/=0)) return
        endif
        if(allocated(final_registers(node)%fine_delta)) then
          if(any(final_registers(node)%fine_delta/=0)) return
        endif
        if(allocated(held_registers(node)%coarse_delta)) then
          if(any(held_registers(node)%coarse_delta/=0)) return
        endif
        if(allocated(held_registers(node)%fine_delta)) then
          if(any(held_registers(node)%fine_delta/=0)) return
        endif
      enddo
      endif
    else if(selected_depth<0.and.prepared_redist) then
      if(tree%topology%level_count()==1) then
        call synchronize_reactive_eb_tree_flux_3d(species,candidate,hf,valid,0,0,boundary_faces,frame)
        if(.not.valid) return
      endif
      do sync_depth=tree%topology%level_count()-2,0,-1
        call synchronize_reactive_eb_tree_register_flux_3d(species,candidate,hf,register_geometry,final_registers, &
          sync_depth,valid,boundary_faces,frame)
        if(.not.valid) return
      enddo
      do node=1,size(tree%nodes)
        if(allocated(final_registers(node)%coarse_delta)) then
          if(any(final_registers(node)%coarse_delta/=0)) return
        endif
        if(allocated(final_registers(node)%fine_delta)) then
          if(any(final_registers(node)%fine_delta/=0)) return
        endif
      enddo
    else
    if(selected_depth<0) then
      if(redistribution_method/='none') before_sync=candidate
      call synchronize_reactive_eb_tree_flux_3d(species,candidate,hf,valid,boundary_faces=boundary_faces,frame=frame)
      if(.not.valid) then
        do sync_depth=tree%topology%level_count()-2,0,-1
          call synchronize_reactive_eb_tree_flux_3d(species,candidate,hf,valid, &
            sync_depth,sync_depth+1,boundary_faces,frame)
          if(.not.valid) then
            call repair_reactive_eb_tree_interface_budget_3d(species,candidate,hf,valid,sync_depth,boundary_faces,frame)
            if(.not.valid) return
            call synchronize_reactive_eb_tree_flux_3d(species,candidate,hf,valid, &
              sync_depth,sync_depth+1,boundary_faces,frame)
          endif
          if(.not.valid) return
        enddo
        call synchronize_reactive_eb_tree_flux_3d(species,candidate,hf,valid,boundary_faces=boundary_faces,frame=frame)
      endif
      if(.not.valid) return
      if(redistribution_method/='none') then
        ! Redistribute the reflux increment too: applying the same linear
        ! operator to rates and flux corrections preserves composite closure.
        sync_rate=hydro
        do node=1,size(tree%nodes)
          sync_rate%fields(node)%values=(candidate%nodes(node)%state-before_sync%nodes(node)%state)/dt
        end do
        if(redistribution_method=='state_redist') then
          call redistribute_eb_tree_state_3d(old,sync_rate,redistributed,valid,redist_target,boundary_faces, &
          max_order=redist_order)
        else
          call redistribute_sdc_eb_tree_flux_rate_3d(old,sync_rate,redistributed,valid,boundary_faces)
        end if
        if(.not.valid) return
        do node=1,size(tree%nodes)
          candidate%nodes(node)%state=before_sync%nodes(node)%state+dt*redistributed%fields(node)%values
        end do
        call recover(valid)
        if(.not.valid) return
        call average_down_reactive_eb_tree_3d(species,candidate,valid,frame)
        if(.not.valid) return
      end if
    end if
    endif
    if(present(internal_energy)) then
      if(selected_depth<0.and..not.present(energy_correction)) then
        if(prepared_redist) then
          if(tree%topology%level_count()==1) then
            call synchronize_reactive_eb_tree_scalar_flux_3d(species,candidate,held_candidate,held_flux,valid, &
              0,0,boundary_faces)
            if(.not.valid) return
          endif
          do sync_depth=tree%topology%level_count()-2,0,-1
            call synchronize_reactive_eb_tree_scalar_register_flux_3d(species,candidate,held_candidate,held_flux, &
              register_geometry,held_registers,sync_depth,valid,boundary_faces)
            if(.not.valid) return
          enddo
          do node=1,size(held_registers)
            if(allocated(held_registers(node)%coarse_delta)) then
              if(any(held_registers(node)%coarse_delta/=0)) return
            endif
            if(allocated(held_registers(node)%fine_delta)) then
              if(any(held_registers(node)%fine_delta/=0)) return
            endif
          enddo
        else
          call synchronize_reactive_eb_tree_scalar_flux_3d(species,candidate,held_candidate,held_flux,valid, &
            boundary_faces=boundary_faces)
          if(.not.valid) return
        endif
      endif
    endif
    if(present(energy_correction)) then
      call recover(valid,3)
      if(.not.valid) return
      if(selected_depth<0.and..not.(mol.and.chemistry.and.redistribution_method=='none')) then
        call average_down_reactive_eb_tree_3d(species,candidate,valid,frame,held_candidate)
        if(.not.valid) return
      endif
    endif
    if(allocated(physical_flux).and.selected_depth>=0) call move_alloc(physical_flux,hf)
    boundary=0
    if(selected_depth<=0) then
    n=tree%topology%nodes(1)%shape
    associate(g=>tree%nodes(1)%geometry,f=>hf(1))
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
    if(.not.all(ieee_is_finite(boundary)).or..not.all(ieee_is_finite(external))) return
    call predictor_field(e1,d1,a1,np)
    if(.not.op%is_valid(tree%topology).or..not.np%is_valid(tree%topology)) return
    if(present(redistribution_transfer)) then
      final_transfer=hydro_transfer
      do node=1,size(tree%nodes)
        if(allocated(admitted_mean_transfer%fields)) then
          final_transfer%fields(node)%values=hydro_transfer%fields(node)%values+ &
            admitted_mean_transfer%fields(node)%values
        else
          final_transfer%fields(node)%values=hydro_transfer%fields(node)%values+ &
            .5_dp*(transfer0%fields(node)%values+transfer1%fields(node)%values)
        endif
      end do
      if(.not.final_transfer%is_valid(tree%topology)) return
    end if
    ! EOS recovery permits roundoff-sized negative trial densities. Published
    ! states must satisfy the strict nonnegative checkpoint contract instead.
    ! Check after reflux/redistribution and before publishing any output.
    do node=1,size(tree%nodes)
      if(.not.active(node)) cycle
      n=tree%topology%nodes(node)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(candidate%nodes(node)%geometry%volume_fraction(i,j,k)<=0) cycle
        if(any(candidate%nodes(node)%state(6:,i,j,k)<0._dp)) return
      end do
      end do
      end do
    end do
    if(present(heat_release).and.selected_depth<0) then
      call average_down_reactive_eb_tree_field_3d(candidate,heat,valid)
      if(.not.valid) return
    endif
    if(present(energy_correction)) then
      numerical_change=initial_reset+final_reset+post_reset+accepted_predictor_reset
      if(.not.ieee_is_finite(energy_correction+numerical_change)) return
    endif
    pressure_work_change=0
    if((present(internal_energy_pressure_work_change).or.godunov_stage_capture).and..not.mol) then
      call field_inventory(held_hydro_work,pressure_work_inventory,valid)
      if(.not.valid) return
      pressure_work_change=-dt*pressure_work_inventory(1)
      if(.not.ieee_is_finite(pressure_work_change)) return
    endif
    if(present(reaction_source)) then
      do node=1,size(tree%nodes)
        if(.not.active(node).or..not.chemistry) cycle
        do c=1,nv-4
          where(old%nodes(node)%geometry%volume_fraction<=0) &
            ir%fields(node)%values(c,:,:,:)=reaction_source%fields(node)%values(c,:,:,:)
        end do
      end do
      if(selected_depth<0) then
        call average_down_reactive_eb_tree_field_3d(candidate,ir,valid)
        if(.not.valid) return
      endif
    end if
    if(present(held_interval_history)) then
      held_history_work%depth=selected_depth;held_history_work%interval=dt
      held_history_work%forcing=rates;held_history_work%mass_source=e0
      held_history_work%local_tree=candidate;held_history_work%local_internal_energy=held_candidate
      held_history_work%reaction_source=ir
      held_history_work%old_flux=f0;held_history_work%spatial_flux=hf
      held_history_work%old_internal_energy_flux=held_flux0;held_history_work%internal_energy_flux=held_flux
      held_history_work%internal_energy_spatial_rate=held_rate0
      do node=1,size(tree%nodes)
        held_history_work%mass_source%fields(node)%values=.5_dp*(e0%fields(node)%values+e1%fields(node)%values)
        if(present(auxiliary_source)) held_history_work%mass_source%fields(node)%values= &
          held_history_work%mass_source%fields(node)%values+.5_dp*(a0%fields(node)%values+a1%fields(node)%values)
        if(allocated(admitted_held_rate%fields)) then
          held_history_work%internal_energy_spatial_rate%fields(node)%values=admitted_held_rate%fields(node)%values
        else
        held_history_work%internal_energy_spatial_rate%fields(node)%values= &
          .5_dp*(held_rate0%fields(node)%values+held_rate1%fields(node)%values)
        endif
        ! Match held_update: accepted MOL funding adds this signed scalar
        ! loan rate to the admitted mean. Its paired faces already contain
        ! the loan; archive both, without relabeling it as a volume source.
        if(mol_source_funded.and.active(node)) then
          n=tree%topology%nodes(node)%shape
          do k=1,n(3)
          do j=1,n(2)
          do i=1,n(1)
            if(all(held_hydro%fields(node)%values(:,i,j,k)==0)) cycle
            call add_eb_tree_delta_cell_3d(held_history_work%internal_energy_spatial_rate% &
              fields(node)%values(:,i,j,k),held_hydro%fields(node)%values(:,i,j,k),valid)
            if(.not.valid) return
          enddo
          enddo
          enddo
        endif
      enddo
      held_history_work%chemical_change=chemical
      held_history_work%wall_change=wall;held_history_work%source_change=external
      allocate(held_history_work%auxiliary_change(naux),source=0._dp)
      if(present(auxiliary_source)) held_history_work%auxiliary_change=auxiliary
      if(held_prepared) held_history_work%prepared_stages=held_stages
      held_history_work%numerical_energy_correction=numerical_change
      held_history_work%forcing_energy_correction=accepted_predictor_reset
      held_history_work%internal_energy_projection_change=projection_inventory(1)
      held_history_work%internal_energy_closure_change=closure_inventory(1)
      held_history_work%initial_internal_energy_change=initial_scalar_reset

    endif
    if(present(source_funding_impulse)) then
      funding_receipt=source_funding_stock
      do node=1,size(tree%nodes)
        funding_receipt%fields(node)%values=0
        if(active(node)) cycle
        funding_receipt%fields(node)%values(:nv,:,:,:)=dt*hydro%fields(node)%values
        funding_receipt%fields(node)%values(nv+1,:,:,:)=dt*held_hydro%fields(node)%values(1,:,:,:)
        n=tree%topology%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(covered(tree%topology,node,[i,j,k])) funding_receipt%fields(node)%values(:,i,j,k)=0
        enddo
        enddo
        enddo
        if(any(.not.ieee_is_finite(funding_receipt%fields(node)%values))) return
      enddo
    endif
    if(present(held_interval_history).and.present(source_funding_stock)) then
      held_history_work%has_source_funding=.true.
      held_history_work%source_funding_stock=source_funding_stock
      held_history_work%source_funding_impulse=funding_receipt
    endif
    if(godunov_stage_capture) then
      call finish_godunov_history(valid)
      if(.not.valid) return
      if(present(godunov_interval_history)) &
        call move_godunov_interval_3d(godunov_history,godunov_interval_history)
    endif
    if(present(held_interval_history)) held_interval_history=held_history_work
    if(present(source_funding_impulse)) call move_alloc(funding_receipt%fields,source_funding_impulse%fields)
    if(present(reaction_source)) reaction_source=ir
    if(present(heat_release)) heat_release=heat
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
    if(present(old_predictor)) old_predictor=op
    if(present(new_predictor)) new_predictor=np
    if(present(auxiliary_change)) auxiliary_change=auxiliary
    if(present(redistribution_transfer)) redistribution_transfer=final_transfer
    tree=candidate;boundary_change=boundary;wall_change=wall;chemical_change=chemical
    source_change=external;minimum_theta=theta;ok=.true.
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
      type(amr_tree_state_3d) :: gas,temp
      type(eb_geometry_3d), allocatable :: geometry(:)
      good=.false.
      if(present(replay_interval).or.native_external_frozen) return
      call history_pack_tree(old,gas,temp,good)
      if(.not.good) return
      allocate(geometry(size(tree%nodes)));geometry=tree%nodes%geometry
      good=valid_native_mol_source_prefix_3d(native_source_replay,tree%topology,geometry,selected_depth,dt,iterations, &
        chemistry,transport_flags,force,[0._dp,0._dp,0._dp],held_policy,gas,temp,held_old,ir, &
        [( .true.,node=1,size(tree%nodes))],energy_correction+initial_reset, &
        present(auxiliary_source),held_auxiliary,frame,source_funding_stock, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected)
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
      if(native_prefix_observed<1.or.godunov_recovery_forbidden.or.godunov_packet_pending) return
      if(.not.allocated(held_history_work%external_stages)) return
      result%depth=selected_depth;result%iterations=iterations
      result%accepted_source_stages=native_prefix_observed;result%interval=dt;result%force=force
      result%separate_momentum=selected_momentum
      result%consistent_velocity=selected_velocity
      result%paired_pressure=selected_pressure
      result%connected_directions=selected_connected
      result%chemistry=chemistry;result%transport_flags=transport_flags
      result%has_auxiliary_source=present(auxiliary_source);result%auxiliary_has_internal_energy=held_auxiliary
      result%has_frame=present(frame);if(present(frame)) result%frame=frame
      result%policy=held_policy;result%topology=tree%topology
      allocate(result%geometry(size(tree%nodes)))
      result%geometry=tree%nodes%geometry
      result%initial_energy_correction=energy_correction+initial_reset
      call history_pack_tree(old,result%initial_gas,result%initial_temperature,ready)
      if(.not.ready) return
      result%initial_internal_energy=held_old
      result%initial_reaction_source=held_history_work%initial_reaction_source
      result%has_source_funding=present(source_funding_stock)
      if(present(source_funding_stock)) result%source_funding_stock=source_funding_stock
      result%external_stages=held_history_work%external_stages
      if(allocated(held_history_work%auxiliary_stages)) result%auxiliary_stages=held_history_work%auxiliary_stages
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
      ready=local_ok;if(.not.ready) return
      local_ok=.not.godunov_recovery_forbidden.and..not.godunov_packet_pending
      ready=local_ok;if(.not.ready) return
      local_ok=godunov_history%accepted_source_stages>0
      if(godunov_history%accepted_source_stages<=iterations) then
        associate(front=>godunov_history%sources(godunov_history%accepted_source_stages))
        local_ok=local_ok.or.front%auxiliary_accepted.or.front%provider_accepted
        end associate
      endif
      ready=local_ok;if(.not.ready) return
      godunov_history%source_prefix_only=.true.
      local_ok=valid_godunov_source_cache_3d(godunov_history,tree%topology,nv,naux,allow_partial=.true.)
      ready=local_ok;if(.not.ready) return
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
      if(present(evaluate_reaction)) then
        call classified_held_prepared_reactor(species,reactions,energy,source_species,source_energy,duration,rtol,atol, &
          density,temperature,steps,ok,maximum_steps)
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
      good=local_ok;if(.not.good) return
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
      good=local_ok;if(.not.good) return
      if(kind==1) then
        call copy_godunov_field_3d(halos,saved%provider_halos,local_ok)
      good=local_ok;if(.not.good) return
      endif
      call copy_godunov_fluxes_3d(faces,saved%provider_faces,local_ok)
      good=local_ok;if(.not.good) return
      saved%provider_kind=kind;saved%provider_limiter=limiter
      saved%provider_accepted=.true.;saved%index=stage
      end associate
      good=.true.
    end subroutine capture_provider_packet
    subroutine validate_godunov_replay(good)
      logical, intent(out) :: good
      logical :: local_ok
      integer :: stage,p,d(3),expected
      good=.false.
      associate(saved=>godunov_replay_interval)
      local_ok=valid_godunov_source_cache_3d(saved,tree%topology,nv,naux,allow_partial=godunov_partial)
      good=local_ok;if(.not.good) return

      local_ok=allocated(saved%topology%nodes).and.allocated(saved%sources)
      good=local_ok;if(.not.good) return
      local_ok=size(saved%topology%nodes)==size(tree%nodes).and.lbound(saved%sources,1)==0.and. &
        ubound(saved%sources,1)==iterations.and.saved%depth==selected_depth.and.saved%iterations==iterations.and. &
        ieee_is_finite(saved%interval).and.all(ieee_is_finite(saved%force))
      good=local_ok;if(.not.good) return
      local_ok=saved%interval==dt.and.all(saved%force==force).and.(saved%chemistry.eqv.chemistry).and. &
        all(saved%transport_flags.eqv.transport_flags).and.(saved%prepared_redistribution.eqv.prepared_redist).and. &
        (saved%has_auxiliary.eqv.present(auxiliary_source)).and. &
        (saved%has_spatial_provider.eqv.present(evaluate_spatial_diffusion)).and. &
        (saved%spatial_diffusion_is_source.eqv.external_diffusion_source).and. &
        (saved%auxiliary_has_internal_energy.eqv.held_auxiliary)
      local_ok=local_ok.and..not.saved%has_valid_diffusion_provider
      good=local_ok;if(.not.good) return
      local_ok=saved%has_source_funding.eqv.present(source_funding_stock)
      good=local_ok;if(.not.good) return
      if(present(source_funding_stock)) then
        local_ok=.true.
        do p=1,size(tree%nodes)
          local_ok=local_ok.and.all(saved%source_funding_stock%fields(p)%values== &
            source_funding_stock%fields(p)%values)
        enddo
        good=local_ok;if(.not.good) return
      endif
      do p=1,size(tree%nodes)
        local_ok=saved%topology%nodes(p)%parent==tree%topology%nodes(p)%parent.and. &
          saved%topology%nodes(p)%level==tree%topology%nodes(p)%level.and. &
          all(saved%topology%nodes(p)%shape==tree%topology%nodes(p)%shape)
      good=local_ok;if(.not.good) return
      enddo
      do stage=0,iterations
        if(.not.(saved%sources(stage)%auxiliary_accepted.or.saved%sources(stage)%provider_accepted).and. &
          stage>=saved%accepted_source_stages) cycle
        local_ok=saved%sources(stage)%index==stage
      good=local_ok;if(.not.good) return
        if(saved%sources(stage)%auxiliary_accepted) then
          call valid_cached_field(saved%sources(stage)%auxiliary,naux,0,good)
          if(.not.good) return
        endif
        expected=0
        if(present(evaluate_spatial_diffusion)) then
          expected=1
          if(external_diffusion_source.and.stage>1.and..not.prepared_redist) expected=0
        endif
        if(.not.saved%sources(stage)%provider_accepted) expected=0
        local_ok=saved%sources(stage)%provider_kind==expected
      good=local_ok;if(.not.good) return
        if(expected==0) cycle
        call valid_cached_field(saved%sources(stage)%provider_rate,nv,0,good)
        if(.not.good) return
        if(expected==1) then
          call valid_cached_field(saved%sources(stage)%provider_halos,nv,6,good)
          if(.not.good) return
        endif
        local_ok=allocated(saved%sources(stage)%provider_faces).and.ieee_is_finite(saved%sources(stage)%provider_limiter)
      good=local_ok;if(.not.good) return
        local_ok=size(saved%sources(stage)%provider_faces)==size(tree%nodes).and. &
          saved%sources(stage)%provider_limiter>=0.and.saved%sources(stage)%provider_limiter<=1
      good=local_ok;if(.not.good) return
        do p=1,size(tree%nodes)
          if(.not.active(p)) cycle
          local_ok=valid_eb_integrated_flux_3d(saved%sources(stage)%provider_faces(p),tree%nodes(p)%geometry,nv)
          if(.not.local_ok) exit
          if(expected==1) then
            d=tree%topology%nodes(p)%shape
            local_ok=all(saved%sources(stage)%provider_halos%fields(p)%values(:,7:d(1)+6,7:d(2)+6,7:d(3)+6)== &
              saved%sources(stage)%provider_rate%fields(p)%values)
            if(.not.local_ok) exit
          endif
        enddo
      good=local_ok;if(.not.good) return
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
      good=local_ok;if(.not.good) return
      local_ok=size(field%fields)==size(tree%nodes).and.lbound(field%fields,1)==1
      good=local_ok;if(.not.good) return
      do p=1,size(tree%nodes)
        if(padding>0.and..not.active(p)) then
          if(allocated(field%fields(p)%values)) local_ok=.false.
          cycle
        endif
        if(.not.allocated(field%fields(p)%values)) then
          local_ok=.false.;exit
        endif
        d=tree%topology%nodes(p)%shape
        local_ok=all(shape(field%fields(p)%values)==[components,d+2*padding]).and. &
          all(lbound(field%fields(p)%values)==1)
        if(.not.local_ok) exit
        local_ok=all(ieee_is_finite(field%fields(p)%values))
        if(.not.local_ok) exit
      enddo
      good=local_ok;if(.not.good) return
      good=.true.
    end subroutine valid_cached_field
    subroutine history_pack_tree(gas,values,temperatures,good)
      type(eb_reactive_tree_3d), intent(in) :: gas
      type(amr_tree_state_3d), intent(out) :: values,temperatures
      logical, intent(out) :: good
      integer :: p,allocation_status,d(3)
      good=.false.
      allocate(values%fields(size(gas%nodes)),temperatures%fields(size(gas%nodes)),stat=allocation_status)
      if(allocation_status/=0) return
      do p=1,size(gas%nodes)
        d=gas%topology%nodes(p)%shape
        allocate(values%fields(p)%values,source=gas%nodes(p)%state,stat=allocation_status)
        if(allocation_status/=0) return
        allocate(temperatures%fields(p)%values(1,d(1),d(2),d(3)),stat=allocation_status)
        if(allocation_status/=0) return
        temperatures%fields(p)%values(1,:,:,:)=gas%nodes(p)%temperature
      enddo
      good=.true.
    end subroutine history_pack_tree
    subroutine begin_godunov_history(good)
      logical, intent(out) :: good
      allocate(godunov_history,stat=history_status)
      good=.false.;if(history_status/=0) return
      allocate(godunov_history%sources(0:iterations),godunov_history%hydro(iterations),stat=history_status)
      if(history_status/=0) return
      godunov_history%depth=selected_depth;godunov_history%interval=dt;godunov_history%iterations=iterations
      godunov_history%force=force;godunov_history%chemistry=chemistry;godunov_history%transport_flags=transport_flags
      godunov_history%prepared_redistribution=prepared_redist;godunov_history%policy=held_policy
      if(present(spatial_flux_weight)) godunov_history%face_weight=spatial_flux_weight
      godunov_history%replayed=present(godunov_replay_interval)
      godunov_history%has_auxiliary=present(auxiliary_source)
      godunov_history%has_source_funding=present(source_funding_stock)
      if(present(source_funding_stock)) godunov_history%source_funding_stock=source_funding_stock
      godunov_history%has_spatial_provider=present(evaluate_spatial_diffusion)
      godunov_history%spatial_diffusion_is_source=external_diffusion_source
      godunov_history%auxiliary_has_internal_energy=held_auxiliary
      if(present(energy_correction)) godunov_history%initial_energy_correction=energy_correction
      godunov_history%topology=tree%topology
      call history_pack_tree(tree,godunov_history%clock_gas,godunov_history%clock_temperature,good)
      if(.not.good) return
      godunov_history%clock_internal_energy=internal_energy
      call history_pack_tree(old,godunov_history%initial_gas,godunov_history%initial_temperature,good)
      if(.not.good) return
      godunov_history%initial_internal_energy=held_old
      godunov_history%initial_reaction_source=ir
      godunov_history%clock_reaction_source=ir
      if(present(ancestor_start)) then
        call history_pack_tree(ancestor_start,godunov_history%ancestor_start, &
          godunov_history%ancestor_temperature_start,good)
        if(.not.good) return
        call history_pack_tree(ancestor_end,godunov_history%ancestor_end,godunov_history%ancestor_temperature_end,good)
        if(.not.good) return
      godunov_history%ancestor_internal_energy_start=ancestor_internal_energy_start
      godunov_history%ancestor_internal_energy_end=ancestor_internal_energy_end
        godunov_history%fractions_start=fractions_start;godunov_history%fractions_end=fractions_end
      endif
      if(present(source_start)) then
      godunov_history%source_start=source_start
      godunov_history%source_end=source_end
      endif
      good=.true.
    end subroutine begin_godunov_history
    subroutine capture_godunov_source(stage,gas,external,diffusion,auxiliary,faces,good)
      integer, intent(in) :: stage
      type(eb_reactive_tree_3d), intent(in) :: gas
      type(amr_tree_state_3d), intent(in) :: external,diffusion,auxiliary
      type(eb_integrated_flux_3d), allocatable, intent(in) :: faces(:)
      logical, intent(out) :: good
      godunov_history%sources(stage)%index=stage
      godunov_history%sources(stage)%coarse_registers_accumulated=prepared_redist.and.stage>1
      associate(saved=>godunov_history%sources(stage))
      call history_pack_tree(gas,saved%gas,saved%temperature,good)
      if(.not.good) return
      saved%external=external
      saved%rate=diffusion
      if(allocated(faces)) saved%faces=faces
      if(present(auxiliary_source)) then
        if(stage==0) then
      saved%auxiliary=auxiliary_packed0
        else
      saved%auxiliary=auxiliary_packed1
        endif
      endif
      if(stage==0) then
      saved%internal_energy=held_old
      saved%scalar_rate=held_rate0
      if(allocated(held_flux0)) saved%scalar_faces=held_flux0
        if(prepared_redist) then
      saved%registers=stage_registers0
      saved%scalar_registers=held_registers0
        endif
      else
      saved%internal_energy=held_candidate
      saved%scalar_rate=held_rate1
      if(allocated(held_flux1)) saved%scalar_faces=held_flux1
        if(prepared_redist) then
      saved%registers=stage_registers1
      saved%scalar_registers=held_registers1
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
        saved%predictor=rates;saved%predictor_rate=hydro
        saved%predictor_scalar_rate=held_hydro;saved%predictor_pressure_work=held_hydro_work
        if(allocated(hf)) saved%predictor_faces=hf
        if(allocated(held_hydro_flux)) saved%predictor_scalar_faces=held_hydro_flux
      endif
      saved%rate=hydro
      saved%scalar_rate=held_hydro
      saved%pressure_work=held_hydro_work
      if(allocated(hf)) saved%faces=hf
      if(allocated(held_hydro_flux)) saved%scalar_faces=held_hydro_flux
      if(prepared_redist) then
      saved%registers=hydro_registers
      saved%scalar_registers=held_hydro_registers
      endif
      end associate
      good=.true.
    end subroutine capture_godunov_hydro
    subroutine finish_godunov_history(good)
      logical, intent(out) :: good
      call history_pack_tree(candidate,godunov_history%final_gas,godunov_history%final_temperature,good)
      if(.not.good) return
      godunov_history%final_internal_energy=held_candidate
      if(present(source_funding_impulse)) godunov_history%source_funding_impulse=funding_receipt
      godunov_history%reaction_source=ir
      if(allocated(hf)) godunov_history%faces=hf
      if(allocated(held_flux)) godunov_history%scalar_faces=held_flux
      if(allocated(registration)) godunov_history%registered_faces=registration
      if(prepared_redist) then
      godunov_history%registers=final_registers
      godunov_history%scalar_registers=held_registers
      endif
      allocate(godunov_history%boundary_change(size(boundary)),godunov_history%wall_change(size(wall)), &
        godunov_history%chemical_change(size(chemical)),godunov_history%source_change(size(external)), &
        godunov_history%auxiliary_change(naux),stat=history_status)
      good=.false.;if(history_status/=0) return
      godunov_history%boundary_change=boundary;godunov_history%wall_change=wall
      godunov_history%chemical_change=chemical;godunov_history%source_change=external
      godunov_history%auxiliary_change=0
      if(present(auxiliary_source)) godunov_history%auxiliary_change=auxiliary
      godunov_history%numerical_energy_change=0
      if(present(energy_correction)) godunov_history%numerical_energy_change=numerical_change
      godunov_history%internal_energy_projection_change=projection_inventory(1)
      godunov_history%internal_energy_closure_change=closure_inventory(1)
      godunov_history%initial_internal_energy_change=initial_scalar_reset
      call copy_godunov_field_3d(op,godunov_history%old_predictor,good)
      if(.not.good) return
      call copy_godunov_field_3d(np,godunov_history%new_predictor,good)
      if(.not.good) return
      godunov_history%pressure_work_change=pressure_work_change
      godunov_history%minimum_theta=theta
      good=.true.
    end subroutine finish_godunov_history
    subroutine redist_boundary_samples(snapshot,good)
      type(eb_reactive_tree_3d), intent(in) :: snapshot
      logical, intent(out) :: good
      real(dp), allocatable :: gt(:,:,:),vf(:,:,:)
      integer :: p,shape3(3)
      good=.false.
      if(allocated(redist_halos)) deallocate(redist_halos)
      if(allocated(redist_geometry)) deallocate(redist_geometry)
      allocate(redist_halos,redist_geometry(size(snapshot%nodes)))
      allocate(redist_halos%fields(size(snapshot%nodes)))
      do p=1,size(snapshot%nodes)
        shape3=snapshot%topology%nodes(p)%shape
        call extract_eb_geometry_3d(extended_geometry(p),[5,5,5],shape3+8,redist_geometry(p),good)
        if(.not.good) return
        allocate(redist_halos%fields(p)%values(nv,shape3(1)+4,shape3(2)+4,shape3(3)+4), &
          gt(shape3(1)+4,shape3(2)+4,shape3(3)+4),vf(shape3(1)+4,shape3(2)+4,shape3(3)+4))
        call fill_reactive_eb_tree_halo_3d(species,snapshot,p,redist_halos%fields(p)%values,gt,good, &
          inflow_state=inflow_state,inflow_mach=inflow_mach,inflow_temperature=inflow_temperature, &
          inflow_pressure=inflow_pressure,boundary_faces=boundary_faces,ghost_volume_fraction=vf, &
          halo_width=2,target_volume_fraction=redist_geometry(p)%volume_fraction, &
          restrict_snapshot=.false.,frame=frame)
        if(.not.good) return
        deallocate(gt,vf)
      enddo
      good=.true.
    end subroutine redist_boundary_samples
    logical function active(p)
      integer, intent(in) :: p
      active=selected_depth<0.or.tree%topology%nodes(p)%level==selected_depth
    end function active
    subroutine validate_held_native_replay(good)
      logical, intent(out) :: good
      good=valid_held_mol_record_3d(tree,replay_interval,nv,frame)
      if(.not.good) return
      if((replay_interval%separate_momentum.neqv.selected_momentum).or. &
          (replay_interval%consistent_velocity.neqv.selected_velocity).or. &
          (replay_interval%paired_pressure.neqv.selected_pressure).or. &
          (replay_interval%connected_directions.neqv.selected_connected)) then
        good=.false.;return
      endif
      good=.not.replay_interval%prepared_redistribution.and.replay_interval%native_iterations==iterations.and. &
        replay_interval%depth==selected_depth.and.replay_interval%interval==dt.and. &
        all(replay_interval%native_force==force).and. &
        (replay_interval%has_auxiliary_source.eqv.present(auxiliary_source)).and. &
        all(replay_interval%native_providers.eqv.[present(evaluate_spatial_diffusion),.false.]).and. &
        replay_interval%native_provider_halo_width==callback_width
      if(.not.good) return
      good=valid_mol_native_provider_cache_3d(tree,replay_interval,nv)
      if(.not.good) return
      good=good.and.replay_interval%policy%eta2==held_policy%eta2.and. &
        (replay_interval%policy%allow_small_energy.eqv.held_policy%allow_small_energy).and. &
        (replay_interval%policy%allow_negative_energy.eqv.held_policy%allow_negative_energy).and. &
        (replay_interval%policy%update_total_from_internal.eqv.held_policy%update_total_from_internal)
    end subroutine validate_held_native_replay

    subroutine validate_held_prepared_replay(good)
      logical, intent(out) :: good
      good=valid_held_prepared_mol_record_3d(tree,replay_interval,nv,selected_depth,dt,iterations, &
        redist_target,redist_order,held_policy,frame,present(auxiliary_source))
    end subroutine validate_held_prepared_replay

    subroutine classified_held_prepared_reactor(species,reactions,initial_energy_density,external_species,external_energy, &
        duration,rtol,atol,species_density,temperature,completed_steps,ok,maximum_steps)
      type(nasa7_species), intent(in) :: species(:)
      type(elementary_reaction), intent(in) :: reactions(:)
      real(dp), intent(in) :: initial_energy_density,external_species(:),external_energy,duration,rtol,atol
      real(dp), intent(inout) :: species_density(:),temperature
      integer, intent(out) :: completed_steps
      logical, intent(out) :: ok
      integer, optional, intent(in) :: maximum_steps
      real(dp) :: probe(size(species_density)),probe_t
      integer :: probe_steps
      logical :: native_ok
      probe=species_density;probe_t=temperature
      call evaluate_reaction(species,reactions,initial_energy_density,external_species,external_energy, &
        duration,rtol,atol,species_density,temperature,completed_steps,ok,maximum_steps)
      if(ok) return
      ! Classification proves admissibility; it never adopts the native result
      ! or retries a supplied provider that refused that admissible input.
      call advance_forced_reactor_interval(species,reactions,initial_energy_density,external_species,external_energy, &
        duration,rtol,atol,probe,probe_t,probe_steps,native_ok,maximum_steps=maximum_steps)
      if(native_ok) prepared_provider_refused=.true.
      godunov_recovery_forbidden=godunov_recovery_forbidden.or.prepared_provider_refused
      if(present(provider_refused)) provider_refused=prepared_provider_refused
    end subroutine classified_held_prepared_reactor

    subroutine run_held_prepared_reaction(good)
      logical, intent(out) :: good
      integer :: p
      held_candidate=held_old
      if(associated(held_prepared_integrator)) then
        call react_sdc_eb_tree_3d(species,reactions,old,rates,dt,rtol,atol,u,t,ir,good,target_depth, &
          evaluate_reaction=held_prepared_integrator,frame=frame,heat_release=heat, &
          internal_energy=held_candidate,old_internal_energy=held_old)
      else
        call react_sdc_eb_tree_3d(species,reactions,old,rates,dt,rtol,atol,u,t,ir,good,target_depth, &
          frame=frame,heat_release=heat,internal_energy=held_candidate,old_internal_energy=held_old)
      endif
      if(.not.good) return
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        candidate%nodes(p)%state=u%fields(p)%values
        candidate%nodes(p)%temperature=t%fields(p)%values(1,:,:,:)
      enddo
      call recover(good,2)
    end subroutine run_held_prepared_reaction

    function held_face_impulse(g,f,x,y,z) result(value)
      type(eb_geometry_3d), intent(in) :: g
      type(eb_integrated_flux_3d), intent(in) :: f
      integer, intent(in) :: x,y,z
      real(dp) :: value(size(f%x,1))
      value=g%dy*g%dz*(g%x_face_fraction(x-1,y,z)*f%x(:,x-1,y,z)-g%x_face_fraction(x,y,z)*f%x(:,x,y,z))+ &
        g%dx*g%dz*(g%y_face_fraction(x,y-1,z)*f%y(:,x,y-1,z)-g%y_face_fraction(x,y,z)*f%y(:,x,y,z))+ &
        g%dx*g%dy*(g%z_face_fraction(x,y,z-1)*f%z(:,x,y,z-1)-g%z_face_fraction(x,y,z)*f%z(:,x,y,z))
    end function held_face_impulse

    subroutine build_held_prepared_low(stage,good)
      integer, intent(in) :: stage
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: ignored,ignored_held,old_halos
      real(dp) :: nonface(nv),volume,a(1),b(1)
      integer :: p,x,y,z
      good=.false.
      if(.not.allocated(held_low_faces)) then
        call compute_mol_eb_tree_hydro_3d(species,old,extended_geometry,dt,1,ignored,held_low_faces,good, &
          boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          target_depth,ancestor_start,ancestor_end,fractions_start,.true.,frame, &
          internal_energy_rate=ignored_held,internal_energy_flux=held_low_scalar_faces,internal_energy_state=held_old, &
          internal_energy_old_state_halos=old_halos,ancestor_internal_energy_start=ancestor_internal_energy_start, &
          ancestor_internal_energy_end=ancestor_internal_energy_end, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected)
        if(.not.good) return
      endif
      held_stages(stage)%low=held_stages(stage)%rate
      held_stages(stage)%scalar_low=held_stages(stage)%scalar_rate
      held_stages(stage)%low_faces=held_low_faces
      held_stages(stage)%scalar_low_faces=held_low_scalar_faces
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        associate(g=>old%nodes(p)%geometry,h=>held_stages(stage))
        do z=1,g%nz
        do y=1,g%ny
        do x=1,g%nx
          h%low%fields(p)%values(:,x,y,z)=0;h%scalar_low%fields(p)%values(:,x,y,z)=0
          if(g%volume_fraction(x,y,z)<=0) cycle
          volume=g%dx*g%dy*g%dz*g%volume_fraction(x,y,z)
          nonface=h%raw%fields(p)%values(:,x+6,y+6,z+6)-held_face_impulse(g,h%faces(p),x,y,z)/volume/dt
          ! Native closed spatial mass/species terms are entirely face based.
          ! Their low family receives no cancellation-created trace inventory.
          nonface(1)=0;nonface(6:)=0
          h%low%fields(p)%values(:,x,y,z)=nonface+held_face_impulse(g,held_low_faces(p),x,y,z)/volume/dt
          a=held_face_impulse(g,h%scalar_faces(p),x,y,z);b=held_face_impulse(g,held_low_scalar_faces(p),x,y,z)
          h%scalar_low%fields(p)%values(1,x,y,z)=h%scalar_raw%fields(p)%values(1,x+6,y+6,z+6)- &
            a(1)/volume/dt+b(1)/volume/dt
        enddo
        enddo
        enddo
        call post_eb_spatial_density_rate_3d(g,h%low%fields(p)%values,good)
        if(.not.good) return
        if(any(.not.ieee_is_finite(h%scalar_low%fields(p)%values))) then
          good=.false.;return
        endif
        end associate
      enddo
      good=.true.
    end subroutine build_held_prepared_low

    subroutine blend_held_prepared_stage(stage,weight,rate,faces,registers,scalar_rate,scalar_faces,scalar_registers,good)
      integer, intent(in) :: stage
      real(dp), intent(in) :: weight
      type(amr_tree_state_3d), intent(inout) :: rate,scalar_rate
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: faces(:),scalar_faces(:)
      type(state_redist_register_3d), allocatable, intent(inout) :: registers(:),scalar_registers(:)
      logical, intent(out) :: good
      integer :: p
      good=.false.
      held_stages(stage)%spatial_weight=weight
      if(weight<1.and..not.allocated(held_stages(stage)%low%fields)) then
        call build_held_prepared_low(stage,good)
        if(.not.good) return
      endif
      rate=held_stages(stage)%rate;scalar_rate=held_stages(stage)%scalar_rate
      faces=held_stages(stage)%faces;scalar_faces=held_stages(stage)%scalar_faces
      registers=held_stages(stage)%registers;scalar_registers=held_stages(stage)%scalar_registers
      if(weight==1) then
        call add_held_auxiliary_rate(scalar_rate,held_stages(stage)%auxiliary_rate,good)
        return
      endif
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        if(weight==0) then
          rate%fields(p)%values=held_stages(stage)%low%fields(p)%values
          scalar_rate%fields(p)%values=held_stages(stage)%scalar_low%fields(p)%values
          faces(p)=held_low_faces(p);scalar_faces(p)=held_low_scalar_faces(p)
        else
          rate%fields(p)%values=(1-weight)*held_stages(stage)%low%fields(p)%values+weight*rate%fields(p)%values
          scalar_rate%fields(p)%values=(1-weight)*held_stages(stage)%scalar_low%fields(p)%values+ &
            weight*scalar_rate%fields(p)%values
          faces(p)%x=(1-weight)*held_low_faces(p)%x+weight*faces(p)%x
          faces(p)%y=(1-weight)*held_low_faces(p)%y+weight*faces(p)%y
          faces(p)%z=(1-weight)*held_low_faces(p)%z+weight*faces(p)%z
          scalar_faces(p)%x=(1-weight)*held_low_scalar_faces(p)%x+weight*scalar_faces(p)%x
          scalar_faces(p)%y=(1-weight)*held_low_scalar_faces(p)%y+weight*scalar_faces(p)%y
          scalar_faces(p)%z=(1-weight)*held_low_scalar_faces(p)%z+weight*scalar_faces(p)%z
        endif
        if(allocated(registers(p)%coarse_delta)) registers(p)%coarse_delta=weight*registers(p)%coarse_delta
        if(allocated(registers(p)%fine_delta)) registers(p)%fine_delta=weight*registers(p)%fine_delta
        if(allocated(scalar_registers(p)%coarse_delta)) &
          scalar_registers(p)%coarse_delta=weight*scalar_registers(p)%coarse_delta
        if(allocated(scalar_registers(p)%fine_delta)) scalar_registers(p)%fine_delta=weight*scalar_registers(p)%fine_delta
        call post_eb_spatial_density_rate_3d(old%nodes(p)%geometry,rate%fields(p)%values,good)
        if(.not.good) return
      enddo
      call add_held_auxiliary_rate(scalar_rate,held_stages(stage)%auxiliary_rate,good)
    end subroutine blend_held_prepared_stage

    subroutine apply_held_prepared_family(last_stage,weight,good)
      integer, intent(in) :: last_stage
      real(dp), intent(in) :: weight
      logical, intent(out) :: good
      call blend_held_prepared_stage(0,weight,d0,f0,stage_registers0,held_rate0,held_flux0,held_registers0,good)
      if(.not.good) return
      call blend_held_prepared_stage(last_stage,weight,d1,f1,stage_registers1, &
        held_rate1,held_flux1,held_registers1,good)
      if(.not.good) return
      held_spatial_weight=weight
      call combined_rates(.false.)
    end subroutine apply_held_prepared_family

    subroutine limit_held_prepared_pass(last_stage,good,reacted)
      integer, intent(in) :: last_stage
      logical, intent(out) :: good,reacted
      type(eb_reactive_tree_3d) :: best
      type(amr_tree_state_3d) :: best_held,best_ir
      type(amr_tree_state_3d), allocatable :: best_heat
      real(dp) :: high(nv),low(nv),fixed(nv),high_mass,low_mass,upper,lower,weight,best_reset
      integer :: pools(nv-5),r,c,first,other,p,x,y,z,n(3),attempt
      logical :: valid,unsafe
      good=.false.;reacted=.false.;unsafe=.false.
      pools=[(c,c=1,nv-5)]
      do r=1,size(reactions)
        if(.not.allocated(reactions(r)%reactant_stoich).or..not.allocated(reactions(r)%product_stoich)) return
        if(size(reactions(r)%reactant_stoich)/=nv-5.or.size(reactions(r)%product_stoich)/=nv-5) return
        first=0
        do c=1,nv-5
          if(reactions(r)%reactant_stoich(c)==reactions(r)%product_stoich(c)) cycle
          if(first==0) then
            first=pools(c)
          else
            other=pools(c)
            where(pools==other) pools=first
          endif
        enddo
      enddo
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        n=tree%topology%nodes(p)%shape
        do z=1,n(3)
        do y=1,n(2)
        do x=1,n(1)
          if(old%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
          high=old%nodes(p)%state(:,x,y,z)+dt*rates%fields(p)%values(:,x,y,z)
          if(any(.not.ieee_is_finite(high))) return
          unsafe=unsafe.or.high(1)<=0
          do c=1,nv-5
            if(pools(c)/=c) cycle
            unsafe=unsafe.or.sum(high(6:),mask=pools==c)<0
          enddo
        enddo
        enddo
        enddo
      enddo
      if(.not.unsafe) then
        call run_held_prepared_reaction(valid)
        if(prepared_provider_refused) return
        if(valid) then
          theta=min(theta,held_spatial_weight);good=.true.;reacted=.true.;return
        endif
      endif
      call build_held_prepared_low(0,valid)
      if(.not.valid) return
      call build_held_prepared_low(last_stage,valid)
      if(.not.valid) return
      upper=held_spatial_weight
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        n=tree%topology%nodes(p)%shape
        do z=1,n(3)
        do y=1,n(2)
        do x=1,n(1)
          if(old%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
          fixed=.5_dp*(e0%fields(p)%values(:,x,y,z)+e1%fields(p)%values(:,x,y,z))
          if(present(auxiliary_source)) fixed=fixed+.5_dp* &
            (a0%fields(p)%values(:,x,y,z)+a1%fields(p)%values(:,x,y,z))
          if(it==1) fixed(5)=fixed(5)+.5_dp*predictor_reset%fields(p)%values(5,x,y,z)/dt
          low=old%nodes(p)%state(:,x,y,z)+dt*(fixed+.5_dp*( &
            held_stages(0)%low%fields(p)%values(:,x,y,z)+held_stages(last_stage)%low%fields(p)%values(:,x,y,z)))
          high=old%nodes(p)%state(:,x,y,z)+dt*(fixed+.5_dp*( &
            held_stages(0)%rate%fields(p)%values(:,x,y,z)+held_stages(last_stage)%rate%fields(p)%values(:,x,y,z)))
          if(any(.not.ieee_is_finite(low)).or.any(.not.ieee_is_finite(high))) return
          if(low(1)<=0) return
          if(high(1)<=0) upper=min(upper,(low(1)/(low(1)-high(1)))*(1-64*epsilon(1._dp)))
          do c=1,nv-5
            if(pools(c)/=c) cycle
            low_mass=sum(low(6:),mask=pools==c);high_mass=sum(high(6:),mask=pools==c)
            if(low_mass<0) return
            if(high_mass<0) upper=min(upper,(low_mass/(low_mass-high_mass))*(1-64*epsilon(1._dp)))
          enddo
        enddo
        enddo
        enddo
      enddo
      ! All gas/scalar stages, faces and CF deltas use one admitted weight.
      ! Rejected chemistry/EOS trials leave no reset or field tail in the result.
      call apply_held_prepared_family(last_stage,0._dp,valid)
      if(.not.valid) return
      call run_held_prepared_reaction(valid)
      if(.not.valid.or.prepared_provider_refused) return
      best=candidate;best_held=held_candidate;best_ir=ir;best_reset=final_reset
      if(present(heat_release)) then
        allocate(best_heat);best_heat=heat
      endif
      lower=0
      do attempt=1,13
        if(attempt==1) then
          weight=upper
        else
          weight=.5_dp*(lower+upper)
        endif
        if(weight==0) exit
        call apply_held_prepared_family(last_stage,weight,valid)
        if(.not.valid) return
        call run_held_prepared_reaction(valid)
        if(prepared_provider_refused) return
        if(valid) then
          lower=weight;best=candidate;best_held=held_candidate;best_ir=ir;best_reset=final_reset
          if(present(heat_release)) best_heat=heat
          if(attempt==1) exit
        else
          upper=weight
        endif
      enddo
      call apply_held_prepared_family(last_stage,lower,valid)
      if(.not.valid) return
      candidate=best;held_candidate=best_held;ir=best_ir;final_reset=best_reset
      if(present(heat_release)) heat=best_heat
      theta=min(theta,lower);good=.true.;reacted=.true.
    end subroutine limit_held_prepared_pass

    subroutine initialize_held_godunov_stage(gas,stage,good,fraction)
      type(eb_reactive_tree_3d), intent(in) :: gas
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
      allocate(faces(size(tree%nodes)))

      local_ok=.true.
      do p=1,size(tree%nodes)

        rate%fields(p)%values=0
        if(.not.active(p)) cycle
        d=tree%topology%nodes(p)%shape
        allocate(faces(p)%x(1,0:d(1),d(2),d(3)),faces(p)%y(1,d(1),0:d(2),d(3)), &
          faces(p)%z(1,d(1),d(2),0:d(3)))

        faces(p)%x=0;faces(p)%y=0;faces(p)%z=0
      enddo
      if(.not.local_ok) return
      call initialize_inactive_faces(faces,1,ready_stage)
      if(.not.ready_stage) return
      if(prepared_godunov) then
        allocate(registers(size(tree%nodes)))

        if(diffusive.or.present(evaluate_spatial_diffusion)) then
                    allocate(old_scalar%fields(size(tree%nodes)),zero_source%fields(size(tree%nodes)))
          do p=1,size(tree%nodes)
            if(.not.active(p)) cycle
            d=tree%topology%nodes(p)%shape
            allocate(ghost(nv,1-6:d(1)+6,1-6:d(2)+6,1-6:d(3)+6),gt(1-6:d(1)+6,1-6:d(2)+6,1-6:d(3)+6), &
              vf(1-6:d(1)+6,1-6:d(2)+6,1-6:d(3)+6),ghost_held(1-6:d(1)+6,1-6:d(2)+6,1-6:d(3)+6))
            call fill_reactive_eb_tree_halo_3d(species,gas,p,ghost,gt,ready_stage, &
              boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
              inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,ghost_volume_fraction=vf, &
              halo_width=6,target_volume_fraction=register_geometry(p)%volume_fraction,restrict_snapshot=.false., &
              ancestor_start=ancestor_start,ancestor_end=ancestor_end,ancestor_fraction=fraction, &
              internal_energy=stored,ghost_internal_energy=ghost_held, &
              ancestor_internal_energy_start=ancestor_internal_energy_start, &
              ancestor_internal_energy_end=ancestor_internal_energy_end)
            if(.not.ready_stage) return
            allocate(old_scalar%fields(p)%values(1,d(1)+12,d(2)+12,d(3)+12), &
              zero_source%fields(p)%values(1,d(1)+12,d(2)+12,d(3)+12),source=0._dp)
            old_scalar%fields(p)%values(1,:,:,:)=ghost_held
            deallocate(ghost,gt,vf,ghost_held)
          enddo
          first=max(0,selected_depth);last=first
          if(selected_depth<0) last=tree%topology%level_count()-1
          do level=first,last
            call compute_state_redist_level_stage_3d(tree%topology,register_geometry,old_scalar,zero_source,1, &
              level,dt,rate,registers,ready_stage,boundary_faces,redist_target,redist_order,register_scale=1._dp)
            if(.not.ready_stage) return
          enddo
        else
          local_ok=.true.
          do p=1,size(tree%nodes)

            if(.not.active(p)) cycle
            call initialize_state_redist_tree_register_3d(tree%topology,p,1,6,registers(p),ready_stage,boundary_faces,1._dp)
            local_ok=local_ok.and.ready_stage
          enddo
          if(.not.local_ok) return
        endif
      endif
      if(stage==0) then
        held_rate0=rate;call move_alloc(faces,held_flux0)
        if(prepared_godunov) call move_alloc(registers,held_registers0)
      else
        if(prepared_godunov) then
          if(stage>1) then
            local_ok=.true.
            do p=1,size(tree%nodes)

              if(.not.active(p)) cycle
              d=tree%topology%nodes(p)%shape
              call accumulate_state_redist_coarse_register_3d(registers(p),held_registers1(p),1,d+12,ready_stage)
              local_ok=local_ok.and.ready_stage
            enddo
            if(.not.local_ok) return
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
      good=.true.
      if(native_mol_weight==1) return
      ! One PCM/Lax-Friedrichs family from the fresh start of this physical
      ! clock. Fresh ancestor endpoints are shared with the native stages.
      if(.not.allocated(held_low_faces)) then
        call compute_mol_eb_tree_hydro_3d(species,old,extended_geometry,dt,1,ignored,held_low_faces,good, &
          boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          target_depth,ancestor_start,ancestor_end,fractions_start,.true.,frame, &
          internal_energy_rate=ignored_held,internal_energy_flux=held_low_scalar_faces,internal_energy_state=held_old, &
          internal_energy_old_state_halos=old_halos, &
          ancestor_internal_energy_start=ancestor_internal_energy_start, &
          ancestor_internal_energy_end=ancestor_internal_energy_end, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected)
        if(.not.good) return
      endif
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        call blend_native_mol_face_pair_3d(old%nodes(p)%geometry,dt,native_mol_weight, &
          held_low_faces(p),held_low_scalar_faces(p),rate%fields(p)%values,faces(p), &
          scalar_rate%fields(p)%values,scalar_faces(p),good, &
          mass_species_from_faces=.not.present(evaluate_spatial_diffusion))
        if(.not.good) return
      enddo
      theta=min(theta,native_mol_weight)
    end subroutine blend_native_mol_stage

    subroutine held_mol_source(stage,snapshot,stored,rate,scalar_flux,scalar_registers,accepted,fraction, &
        gas_rate,gas_flux,gas_spatial_halos,gas_old_halos,scalar_raw_halos,scalar_old_halos)
      integer, intent(in) :: stage
      type(eb_mol_diagnostics_3d), allocatable :: captured(:)
      real(dp) :: clock_start
      type(eb_reactive_tree_3d), intent(in) :: snapshot
      type(amr_tree_state_3d), intent(in) :: stored
      type(amr_tree_state_3d), intent(inout) :: rate
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: scalar_flux(:)
      type(state_redist_register_3d), allocatable, intent(inout) :: scalar_registers(:)
      logical, intent(out) :: accepted
      real(dp), optional, intent(in) :: fraction(:)
      type(amr_tree_state_3d), intent(inout) :: gas_rate
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: gas_flux(:)
      type(amr_tree_state_3d), optional, intent(inout) :: gas_spatial_halos,gas_old_halos
      type(amr_tree_state_3d), optional, intent(inout) :: scalar_raw_halos
      type(amr_tree_state_3d), optional, intent(inout) :: scalar_old_halos
      type(amr_tree_state_3d) :: old_halos,spatial_halos
      integer :: depth,first_depth,last_depth
      if(present(observe_mol)) then
      call compute_mol_eb_tree_hydro_3d(species,snapshot,extended_geometry,dt,2,gas_rate,gas_flux,accepted, &
        boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        lax_friedrichs=mol_lax_friedrichs,frame=frame,internal_energy_rate=rate, &
        internal_energy_flux=scalar_flux,internal_energy_state=stored,internal_energy_old_state_halos=old_halos, &
        internal_energy_spatial_rate_halos=spatial_halos, &
        spatial_rate_halos=gas_spatial_halos,old_state_halos=gas_old_halos, &
        target_depth=target_depth,ancestor_start=ancestor_start,ancestor_end=ancestor_end, &
        ancestor_fraction=fraction, &
        ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
        diagnostics=captured)
      else
      call compute_mol_eb_tree_hydro_3d(species,snapshot,extended_geometry,dt,2,gas_rate,gas_flux,accepted, &
        boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        lax_friedrichs=mol_lax_friedrichs,frame=frame,internal_energy_rate=rate, &
        internal_energy_flux=scalar_flux,internal_energy_state=stored,internal_energy_old_state_halos=old_halos, &
        internal_energy_spatial_rate_halos=spatial_halos, &
        spatial_rate_halos=gas_spatial_halos,old_state_halos=gas_old_halos, &
        target_depth=target_depth,ancestor_start=ancestor_start,ancestor_end=ancestor_end, &
        ancestor_fraction=fraction, &
        ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected)
      endif
      if(.not.accepted) return
      if(present(observe_mol)) then
        clock_start=0
        if(present(diagnostic_time)) clock_start=diagnostic_time
        call observe_mol(stage,selected_depth,clock_start,dt,snapshot,captured,gas_flux)
      endif
      if(present(scalar_raw_halos)) scalar_raw_halos=spatial_halos
      if(present(scalar_old_halos)) scalar_old_halos=old_halos
      call initialize_inactive_faces(scalar_flux,1,accepted)
      if(.not.accepted) return
      if(.not.prepared_redist) return
      first_depth=max(0,selected_depth);last_depth=first_depth
      if(selected_depth<0) last_depth=tree%topology%level_count()-1
      do depth=first_depth,last_depth
        call compute_state_redist_level_stage_3d(snapshot%topology,extended_geometry,old_halos,spatial_halos,1, &
          depth,dt,rate,scalar_registers,accepted,boundary_faces,redist_target,redist_order,register_scale=.5_dp)
        if(.not.accepted) return
      enddo
    end subroutine held_mol_source
    ! Fund a represented forward source debt from the actual old material
    ! family. All gas components and independent UEINT share one packet;
    ! numerical face transfers and the native hydro RHS change together.
    ! Selected clocks and prepared registers retain their existing replay path.
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
      do p=1,size(tree%nodes)
        matches=matches.and.all(old%nodes(p)%state==saved%initial_gas%fields(p)%values)
        matches=matches.and.all(old%nodes(p)%temperature== &
          saved%initial_temperature%fields(p)%values(1,:,:,:))
        matches=matches.and.all(held_old%fields(p)%values==saved%initial_internal_energy%fields(p)%values)
      enddo
      if(.not.matches) return
      hydro=saved%hydro(it)%rate;held_hydro=saved%hydro(it)%scalar_rate
      held_hydro_work=saved%hydro(it)%pressure_work
      hf=saved%hydro(it)%faces;held_hydro_flux=saved%hydro(it)%scalar_faces
      source_funded=saved%hydro(it)%source_funded
      if(godunov_stage_capture) call capture_godunov_hydro(it,good,retain_predictor=.true.)
      end associate
    end subroutine restore_cached_funding_hydro

    subroutine admit_held_native_corrector(good)
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: trial_h
      real(dp) :: post(reactive_nvar(size(species))),qp(reactive_nprim(size(species))),tg,sound,bias
      integer :: p,d(3),x,y,z
      logical :: unsafe,cell_ok
      good=.false.;unsafe=.false.
      call restore_cached_funding_hydro(good)
      if(.not.good) return
      call held_update(held_old,held_rate0,held_rate1,.5_dp,trial_h,cell_ok)
      if(.not.cell_ok) return
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        d=tree%topology%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(old%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
          ! Virtual coarse sources own no physical inventory. Their reacting
          ! forcing is reconstructed from the admitted fine stage below.
          if(chemistry.and.covered(tree%topology,p,[x,y,z])) cycle
          post=old%nodes(p)%state(:,x,y,z)+ &
            dt*(rates%fields(p)%values(:,x,y,z)+hydro%fields(p)%values(:,x,y,z))
          cell_ok=all(ieee_is_finite(post)).and.post(1)>0
          if(.not.chemistry) cell_ok=cell_ok.and.all(post(6:)>=0)
          ! Signed species forcing can become physical through reaction.
          ! A caloric preview requires positive composition, while the private
          ! built-in trial retains signed forcing. Chemistry cannot fix rho.
          if(cell_ok.and.chemistry.and..not.present(evaluate_reaction)) then
            bias=0
            call compute_reactive_eb_reactor_bias_3d(old%nodes(p)%geometry,[x,y,z],old%nodes(p)%state(:,x,y,z), &
              held_old%fields(p)%values(1,x,y,z),bias,cell_ok,frame)
            if(cell_ok) then
              if(all(post(6:)>=0)) then
              call admit_reactive_eb_reactor_energy_3d(species,old%nodes(p)%geometry,[x,y,z],post,bias, &
                old%nodes(p)%temperature(x,y,z),cell_ok,frame)
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
          ! Previously funded material retains its independently conserved H.
          ! Ordinary unfunded thermal resets keep the pinned policy path.
          if(cell_ok.and.source_funded.and..not.chemistry) then
            call recover_reactive_eb_cell_3d(species,old%nodes(p)%geometry,[x,y,z],post, &
              old%nodes(p)%temperature(x,y,z),qp,tg,sound,cell_ok,frame,trial_h%fields(p)%values(1,x,y,z))
          endif
          unsafe=unsafe.or..not.cell_ok
        enddo
        enddo
        enddo
      enddo
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
      do p=1,size(tree%nodes)
        endpoint%fields(p)%values=old%nodes(p)%state+dt*rates%fields(p)%values
        d=tree%topology%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(.not.covered(tree%topology,p,[x,y,z])) cycle
          ! A restored funded hydro family has no mutable funding flag. A
          ! deficient virtual density still selects the same reconstruction.
          reconstruct=reconstruct.or.endpoint%fields(p)%values(1,x,y,z)<=0
        enddo
        enddo
        enddo
      enddo
      good=.true.
      if(.not.reconstruct) return
      ! Restrict the forced gas vector, not a caloric preview: reaction may
      ! make its composition admissible. The actual old UEINT clock remains
      ! the reactor's authority, and physical leaf rates are never rewritten.
      call average_down_reactive_eb_tree_field_3d(old,endpoint,good)
      if(.not.good) return
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        d=tree%topology%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(.not.covered(tree%topology,p,[x,y,z])) cycle
          rates%fields(p)%values(:,x,y,z)=(endpoint%fields(p)%values(:,x,y,z)- &
            old%nodes(p)%state(:,x,y,z))/dt
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
      initial(:,1,1,1)=old%nodes(node)%state(:,cell(1),cell(2),cell(3))
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
      call compute_reactive_eb_reactor_bias_3d(old%nodes(node)%geometry,cell,initial(:,1,1,1),0._dp, &
        old_minus_projection,finite_projection,frame)
      if(.not.finite_projection) return
      call compute_reactive_eb_reactor_bias_3d(old%nodes(node)%geometry,cell,trial(:,1,1,1),0._dp, &
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
      initial_t=old%nodes(node)%temperature(cell(1),cell(2),cell(3));trial_t=initial_t
      old_h=held_old%fields(node)%values(1,cell(1),cell(2),cell(3));trial_h=old_h
      if(present(frame)) then
        allocate(potential(1,1,1))
        associate(g=>old%nodes(node)%geometry)
        position=[g%x_lower,g%y_lower,g%z_lower]+(real(cell,dp)-.5_dp)*[g%dx,g%dy,g%dz]
        end associate
        call rotating_potential_3d(frame,position,potential(1,1,1),good)
        if(.not.good) return
      endif
      call react_sdc_state_3d(species,reactions,initial,initial_t,forcing,dt,rtol,atol,trial,trial_t,source,good, &
        rotating_potential=potential,internal_energy_density=trial_h,old_internal_energy_density=old_h)
    end subroutine try_default_reactor_cell

    ! Native selected hydro leaves inactive scalar faces unallocated.
    ! A donor transaction can touch those clocks, so represent their absent
    ! native operator with zero faces before the paired inventory/flux check.
    subroutine ensure_selected_scalar_faces(good)
      logical, intent(out) :: good
      call initialize_inactive_faces(hf,nv,good)
      if(.not.good) return
      call initialize_inactive_faces(held_hydro_flux,1,good)
      if(.not.good) return
      ! The local mean-flux limiter allocates selected clocks only. Its
      ! admitted family replaces held_flux0, including otherwise absent
      ! inactive donors, before the complete paired inventory check.
      if(allocated(admitted_held_faces)) &
        call initialize_inactive_faces(admitted_held_faces,1,good)
    end subroutine ensure_selected_scalar_faces

    subroutine initialize_inactive_faces(faces,components,good)
      type(eb_integrated_flux_3d), intent(inout) :: faces(:)
      integer, intent(in) :: components
      logical, intent(out) :: good
      integer :: p,d(3),allocation_status
      good=.true.
      if(selected_depth<0.or..not.native_held_funding) return
      do p=1,size(tree%nodes)
        if(active(p)) cycle
        d=tree%topology%nodes(p)%shape
        if(.not.allocated(faces(p)%x)) then
          allocate(faces(p)%x(components,0:d(1),d(2),d(3)),source=0._dp,stat=allocation_status)
          if(allocation_status/=0) then
            good=.false.;return
          endif
        endif
        if(.not.allocated(faces(p)%y)) then
          allocate(faces(p)%y(components,d(1),0:d(2),d(3)),source=0._dp,stat=allocation_status)
          if(allocation_status/=0) then
            good=.false.;return
          endif
        endif
        if(.not.allocated(faces(p)%z)) then
          allocate(faces(p)%z(components,d(1),d(2),0:d(3)),source=0._dp,stat=allocation_status)
          if(allocation_status/=0) then
            good=.false.;return
          endif
        endif
      enddo
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
      if(require_post.and..not.covered(tree%topology,node,cell)) then
        ! Preserve the native rate/impulse association exactly. A packet
        ! accepted by pre+debt can fail composition closure after cancellation
        ! in this expression, so validate it while the proposal is private.
        next_rate=hydro%fields(node)%values(:,x,y,z)+(pre(:nv)-old%nodes(node)%state(:,x,y,z))/dt
        next_h=held_hydro%fields(node)%values(1,x,y,z)+ &
          (pre(nv+1)-held_old%fields(node)%values(1,x,y,z))/dt
        cell_ok=all(ieee_is_finite(next_rate)).and.ieee_is_finite(next_h)
        if(cell_ok) then
          delivered=old%nodes(node)%state(:,x,y,z)+dt*(rates%fields(node)%values(:,x,y,z)+next_rate)
          stored=funding_raw_h%fields(node)%values(1,x,y,z)+ &
            dt*(next_h-held_hydro%fields(node)%values(1,x,y,z))
          if(present(source_funding_stock).and..not.active(node)) then
            delivered=delivered+source_funding_stock%fields(node)%values(:nv,x,y,z)
            stored=stored+source_funding_stock%fields(node)%values(nv+1,x,y,z)
          endif
          call recover_reactive_eb_cell_3d(species,old%nodes(node)%geometry,cell,delivered, &
            old%nodes(node)%temperature(x,y,z),pv,temp,sound,cell_ok,frame,stored)
          if(cell_ok.and.chemistry) then
            call compute_reactive_eb_reactor_bias_3d(old%nodes(node)%geometry,cell,old%nodes(node)%state(:,x,y,z), &
              held_old%fields(node)%values(1,x,y,z),bias,cell_ok,frame)
            if(cell_ok) call admit_reactive_eb_reactor_energy_3d(species,old%nodes(node)%geometry,cell,delivered, &
              bias,old%nodes(node)%temperature(x,y,z),cell_ok,frame)
            ! A static gas/H endpoint can leave an active donor unable to
            ! complete its actual forced reaction. Admit the same native rate
            ! with a private built-in trial for the accepted corrector family.
            ! A predictor supplies only an unreacted spatial/source estimate.
            ! Inactive clocks carry deferred receipts; caller-owned reactors
            ! retain their ordinary single attempt and are never previewed.
            if(cell_ok.and.funding_corrector_stage.and.active(node).and..not.present(evaluate_reaction)) &
              call try_default_reactor_cell(node,cell,cell_ok,rates%fields(node)%values(:,x,y,z)+next_rate)
          endif
        endif
      endif
      good=cell_ok
    end subroutine admit_native_funding_cell

    subroutine fund_held_native_stage(good,corrector)
      logical, intent(out) :: good
      logical, optional, intent(in) :: corrector
      type(amr_tree_state_3d), allocatable :: reactor_bias
      logical :: final_stage
      type(eb_reactive_tree_3d) :: funded
      type(amr_tree_state_3d) :: gas_debt,held_debt,funded_h,raw_h,next_hydro,next_held
      type(eb_integrated_flux_3d), allocatable :: gas_faces(:),scalar_faces(:),gas_stage_faces(:),held_stage_faces(:)
      real(dp) :: post(reactive_nvar(size(species))),stored,pv(reactive_nprim(size(species))),temp,sound
      integer :: p,d(3),x,y,z,allocation_status
      logical :: cell_ok
      final_stage=.false.
      if(present(corrector)) final_stage=corrector
      call ensure_selected_scalar_faces(good)
      if(.not.good) return
      good=.false.;funded=old;funded_h=held_old;gas_faces=hf;scalar_faces=held_hydro_flux
      if(mol.or.diffusive) then
        gas_stage_faces=f0;held_stage_faces=held_flux0
        if(final_stage) then
          do p=1,size(tree%nodes)
            gas_stage_faces(p)%x=.5_dp*(f0(p)%x+f1(p)%x)
            gas_stage_faces(p)%y=.5_dp*(f0(p)%y+f1(p)%y)
            gas_stage_faces(p)%z=.5_dp*(f0(p)%z+f1(p)%z)
            held_stage_faces(p)%x=.5_dp*(held_flux0(p)%x+held_flux1(p)%x)
            held_stage_faces(p)%y=.5_dp*(held_flux0(p)%y+held_flux1(p)%y)
            held_stage_faces(p)%z=.5_dp*(held_flux0(p)%z+held_flux1(p)%z)
          enddo
        endif
        do p=1,size(tree%nodes)
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
        do p=1,size(tree%nodes)
          if(active(p)) cycle
          raw_h%fields(p)%values=held_old%fields(p)%values+dt*held_hydro%fields(p)%values
        enddo
      endif
      allocate(gas_debt%fields(size(tree%nodes)),held_debt%fields(size(tree%nodes)))
      do p=1,size(tree%nodes)
        gas_debt%fields(p)%values=dt*(rates%fields(p)%values+hydro%fields(p)%values)
        held_debt%fields(p)%values=raw_h%fields(p)%values-held_old%fields(p)%values
        if(present(source_funding_stock).and..not.active(p)) then
          gas_debt%fields(p)%values=gas_debt%fields(p)%values+source_funding_stock%fields(p)%values(:nv,:,:,:)
          held_debt%fields(p)%values(1,:,:,:)=held_debt%fields(p)%values(1,:,:,:)+ &
            source_funding_stock%fields(p)%values(nv+1,:,:,:)
        endif
        ! Covered operator values are virtual; physical fine sources already
        ! carry their own debt. Re-reflux must not charge those values twice.
        d=tree%topology%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(.not.covered(tree%topology,p,[x,y,z])) cycle
          gas_debt%fields(p)%values(:,x,y,z)=0
          held_debt%fields(p)%values(:,x,y,z)=0
        enddo
        enddo
        enddo
      enddo
      if(chemistry) then
        allocate(reactor_bias,stat=allocation_status)
        if(allocation_status/=0) return
        reactor_bias=held_old
        do p=1,size(tree%nodes)
          reactor_bias%fields(p)%values=0
          d=tree%topology%nodes(p)%shape
          do z=1,d(3)
          do y=1,d(2)
          do x=1,d(1)
            if(old%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
            call compute_reactive_eb_reactor_bias_3d(old%nodes(p)%geometry,[x,y,z],old%nodes(p)%state(:,x,y,z), &
              held_old%fields(p)%values(1,x,y,z),reactor_bias%fields(p)%values(1,x,y,z),cell_ok,frame)
            if(.not.cell_ok) return
          enddo
          enddo
          enddo
        enddo

      endif
      funding_raw_h=raw_h
      funding_corrector_stage=final_stage
      call repair_reactive_eb_tree_held_interface_budget_3d(species,funded,funded_h,gas_faces,scalar_faces, &
        cell_ok,0,boundary_faces,frame,gas_debt,held_debt,reactor_energy_bias=reactor_bias, &
        native_projection=admit_native_funding_cell)
      if(.not.cell_ok) return
      next_hydro=hydro;next_held=held_hydro
      do p=1,size(tree%nodes)
        next_hydro%fields(p)%values=hydro%fields(p)%values+ &
          (funded%nodes(p)%state-old%nodes(p)%state)/dt
        next_held%fields(p)%values=held_hydro%fields(p)%values+ &
          (funded_h%fields(p)%values-held_old%fields(p)%values)/dt
        if(.not.all(ieee_is_finite(next_hydro%fields(p)%values)).or. &
            .not.all(ieee_is_finite(next_held%fields(p)%values))) return
        d=tree%topology%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(old%nodes(p)%geometry%volume_fraction(x,y,z)<=0.or. &
            covered(tree%topology,p,[x,y,z])) cycle
          post=old%nodes(p)%state(:,x,y,z)+dt*(rates%fields(p)%values(:,x,y,z)+ &
            next_hydro%fields(p)%values(:,x,y,z))
          stored=raw_h%fields(p)%values(1,x,y,z)+ &
            dt*(next_held%fields(p)%values(1,x,y,z)-held_hydro%fields(p)%values(1,x,y,z))
          if(present(source_funding_stock).and..not.active(p)) then
            post=post+source_funding_stock%fields(p)%values(:nv,x,y,z)
            stored=stored+source_funding_stock%fields(p)%values(nv+1,x,y,z)
          endif
          call recover_reactive_eb_cell_3d(species,old%nodes(p)%geometry,[x,y,z],post, &
            old%nodes(p)%temperature(x,y,z),pv,temp,sound,cell_ok,frame,stored)
          if(.not.cell_ok) return
          if(chemistry) then
            call admit_reactive_eb_reactor_energy_3d(species,old%nodes(p)%geometry,[x,y,z],post, &
              reactor_bias%fields(p)%values(1,x,y,z),old%nodes(p)%temperature(x,y,z),cell_ok,frame)
            if(.not.cell_ok) return
          endif
        enddo
        enddo
        enddo
      enddo
      if(mol.or.diffusive) then
        do p=1,size(tree%nodes)
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
      hydro=next_hydro;held_hydro=next_held;hf=gas_faces;held_hydro_flux=scalar_faces;good=.true.
    end subroutine fund_held_native_stage

    subroutine measure_held_projection(good)
      logical, intent(out) :: good
      integer :: p,x,y,z,d(3)
      logical :: finite,cell_ok
      call held_update(held_old,held_rate0,held_rate1,.5_dp,projection_transport,good)
      if(.not.good) return
      projection_delta=held_old;finite=.true.
      do p=1,size(tree%nodes)
        projection_delta%fields(p)%values=0
        if(.not.active(p)) cycle
        d=tree%topology%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(old%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
          if(covered(tree%topology,p,[x,y,z])) cycle
          call compute_reactive_eb_held_projection(internal_energy%fields(p)%values(1,x,y,z), &
            held_old%fields(p)%values(1,x,y,z),projection_transport%fields(p)%values(1,x,y,z), &
            held_candidate%fields(p)%values(1,x,y,z),projection_delta%fields(p)%values(1,x,y,z),cell_ok)
          finite=finite.and.cell_ok
        enddo
        enddo
        enddo
      enddo
      good=finite
      if(.not.good) return
      call field_inventory(projection_delta,projection_inventory,good)
    end subroutine measure_held_projection

    subroutine held_update(initial,first,second,weight,output,accepted)
      type(amr_tree_state_3d), intent(in) :: initial,first,second
      real(dp), intent(in) :: weight
      type(amr_tree_state_3d), intent(inout) :: output
      logical, intent(out) :: accepted
      integer :: p
      output=initial;accepted=.true.
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        if(allocated(admitted_held_rate%fields)) then
          output%fields(p)%values=initial%fields(p)%values+dt*admitted_held_rate%fields(p)%values
        else
        output%fields(p)%values=initial%fields(p)%values+weight*dt* &
          (first%fields(p)%values+second%fields(p)%values)
        endif
        if(.not.mol.or.mol_source_funded) &
          output%fields(p)%values=output%fields(p)%values+dt*held_hydro%fields(p)%values
        if(.not.all(ieee_is_finite(output%fields(p)%values))) accepted=.false.
      enddo
    end subroutine held_update
    subroutine predictor_field(ext,diff,aux,output)
      type(amr_tree_state_3d), intent(in) :: ext,diff,aux
      type(amr_tree_state_3d), intent(out) :: output
      integer :: p
      output=ext
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        output%fields(p)%values=ext%fields(p)%values+diff%fields(p)%values
        if(present(auxiliary_source)) output%fields(p)%values=output%fields(p)%values+aux%fields(p)%values
        output%fields(p)%values(6:,:,:,:)=output%fields(p)%values(6:,:,:,:)+ir%fields(p)%values(:nv-5,:,:,:)
        output%fields(p)%values(5,:,:,:)=output%fields(p)%values(5,:,:,:)+ir%fields(p)%values(nv-4,:,:,:)
      end do
    end subroutine predictor_field
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
      integer :: p,d(3),x,y,z
      logical :: safe,local_safe,periodic(3)
      good=.false.;local_safe=.true.
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        d=tree%topology%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(old%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
          if(source_funded.and.covered(tree%topology,p,[x,y,z])) cycle
          if(mol_source_funded) then
            ! Funding admits this represented native state. A separate
            ! addition of opposing rates can manufacture a negative trace
            ! and discard an already admissible conservative loan family.
            endpoint=old%nodes(p)%state(:,x,y,z)+dt* &
              (rates%fields(p)%values(:,x,y,z)+hydro%fields(p)%values(:,x,y,z))
          else
            endpoint=old%nodes(p)%state(:,x,y,z)+dt*rates%fields(p)%values(:,x,y,z)
          endif
          if(redistribution_method/='flux_redist') &
            endpoint(6:)=endpoint(6:)+dt*ir%fields(p)%values(:nv-5,x,y,z)
          local_safe=local_safe.and.all(ieee_is_finite(endpoint)).and.endpoint(1)>0.and.all(endpoint(6:)>=0)
        enddo
        enddo
        enddo
      enddo
      if(redistribution_method=='flux_redist'.and.present(redistribution_transfer)) then
        do p=1,size(tree%nodes)
          if(active(p)) cycle
          d=tree%topology%nodes(p)%shape
          do z=1,d(3)
          do y=1,d(2)
          do x=1,d(1)
            if(old%nodes(p)%geometry%volume_fraction(x,y,z)<=0.or.covered(tree%topology,p,[x,y,z])) cycle
            endpoint=old%nodes(p)%state(:,x,y,z)+ &
              .5_dp*(transfer0%fields(p)%values(:,x,y,z)+transfer1%fields(p)%values(:,x,y,z))
            local_safe=local_safe.and.all(ieee_is_finite(endpoint)).and.endpoint(1)>0.and.all(endpoint(6:)>=0)
          enddo
          enddo
          enddo
        enddo
      endif
      safe=local_safe
      if(safe) then
        good=.true.;return
      endif
      ! A nonunit root trial owns one family for the whole forest. Paired
      ! funding has already run; do not hide a second local flux blend here.
      if(native_mol_weight/=1) return
      if(redistribution_method=='flux_redist') then
        call limit_mol_flux_redist_pass(good)
        return
      endif
      fixed_rate=e0;new_rate=d1
      if(present(internal_energy)) then
      call compute_mol_eb_tree_hydro_3d(species,old,extended_geometry,dt,1,low_rate,baseline,safe, &
        boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        target_depth,ancestor_start,ancestor_end,fractions_start,mol_lax_friedrichs,frame, &
        internal_energy_rate=low_held_rate,internal_energy_flux=low_held_faces,internal_energy_state=held_old, &
        internal_energy_old_state_halos=old_held_halos, &
        ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end, &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected)
      else
      call compute_mol_eb_tree_hydro_3d(species,old,mol_geometry,dt,1,low_rate,baseline,safe, &
        boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        target_depth,ancestor_start,ancestor_end,fractions_start,mol_lax_friedrichs,frame)
      endif
      if(.not.safe) return
      allocate(new_flux(size(tree%nodes)))
      if(present(internal_energy)) then
        new_held_rate=held_rate1
        allocate(new_held_faces(size(tree%nodes)))
        do p=1,size(tree%nodes)
          if(.not.active(p)) cycle
          new_held_rate%fields(p)%values=.5_dp*(held_rate0%fields(p)%values+held_rate1%fields(p)%values)
          new_held_faces(p)=held_flux1(p)
          new_held_faces(p)%x=.5_dp*(held_flux0(p)%x+held_flux1(p)%x)
          new_held_faces(p)%y=.5_dp*(held_flux0(p)%y+held_flux1(p)%y)
          new_held_faces(p)%z=.5_dp*(held_flux0(p)%z+held_flux1(p)%z)
        enddo
      endif
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
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
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        d=tree%topology%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(old%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
          endpoint=old%nodes(p)%state(:,x,y,z)+dt*(low_rate%fields(p)%values(:,x,y,z)+fixed_rate%fields(p)%values(:,x,y,z))
          local_safe=local_safe.and.all(ieee_is_finite(endpoint)).and.endpoint(1)>0.and.all(endpoint(6:)>=0)
        enddo
        enddo
        enddo
      enddo
      safe=local_safe
      if(.not.safe) then
        do p=1,size(tree%nodes)
          if(.not.active(p)) cycle
          fixed_rate%fields(p)%values=.5_dp*(e0%fields(p)%values+e1%fields(p)%values)
          if(present(auxiliary_source)) fixed_rate%fields(p)%values=fixed_rate%fields(p)%values+ &
            .5_dp*(a0%fields(p)%values+a1%fields(p)%values)
          if(mol_source_funded) fixed_rate%fields(p)%values=fixed_rate%fields(p)%values+hydro%fields(p)%values
        enddo
      endif
      periodic=[boundary_faces(1)=='periodic',boundary_faces(3)=='periodic',boundary_faces(5)=='periodic']
      if(present(internal_energy)) then
      call apply_eb_tree_flux_correction_3d(old,fixed_rate,baseline,dt,new_rate,new_flux, &
        safe,periodic,target_depth=selected_depth, &
        internal_energy_baseline=low_held_faces,internal_energy_rate=new_held_rate,internal_energy_flux=new_held_faces)
        if(.not.safe) return
        call move_alloc(new_held_faces,admitted_held_faces)
        admitted_held_rate=new_held_rate
      else
      call apply_eb_tree_flux_correction_3d(old,fixed_rate,baseline,dt,new_rate,new_flux, &
        safe,periodic,target_depth=selected_depth)
      endif
      if(.not.safe) return
      admitted_mean_rate=new_rate;admitted_mean_faces=new_flux
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
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
      residual=source;good=.true.
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        associate(g=>old%nodes(p)%geometry,f=>faces(p))
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
            if(c==1.or.c>=6) then
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

    subroutine limit_mol_flux_redist_pass(good)
      logical, intent(out) :: good
      type(amr_tree_state_3d) :: fixed_rate,low_rate,low_transfer,raw_rate
      type(eb_integrated_flux_3d), allocatable :: low_flux(:)
      real(dp) :: low(nv),high(nv),blend,weight,divergence(nv),volume
      integer :: p,d(3),x,y,z,c,retreat
      logical :: valid,safe
      good=.false.;fixed_rate=e0
      call compute_mol_eb_tree_hydro_3d(species,old,mol_geometry,dt,1,raw_rate,low_flux,safe, &
        boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        target_depth,ancestor_start,ancestor_end,fractions_start,.true.,frame)
      if(.not.safe) return
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        associate(g=>old%nodes(p)%geometry,f=>low_flux(p))
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
        do p=1,size(tree%nodes)
          low_transfer%fields(p)%values=0
        enddo
      endif
      valid=.true.
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        call post_eb_spatial_density_rate_3d(old%nodes(p)%geometry,low_rate%fields(p)%values,valid)
        if(.not.valid) exit
      enddo
      safe=valid
      if(.not.safe) return
      ! Only accepted external mass can fund transport. I_R predicts a
      ! future chemical increment and is not an available species inventory.
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        fixed_rate%fields(p)%values=.5_dp*(e0%fields(p)%values+e1%fields(p)%values)
        if(present(auxiliary_source)) fixed_rate%fields(p)%values=fixed_rate%fields(p)%values+ &
          .5_dp*(a0%fields(p)%values+a1%fields(p)%values)
      enddo
      valid=.true.;blend=1
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        d=tree%topology%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(old%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
          low=old%nodes(p)%state(:,x,y,z)+ &
            dt*(fixed_rate%fields(p)%values(:,x,y,z)+low_rate%fields(p)%values(:,x,y,z))
          high=old%nodes(p)%state(:,x,y,z)+dt*(fixed_rate%fields(p)%values(:,x,y,z)+ &
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
        do p=1,size(tree%nodes)
          if(active(p)) cycle
          d=tree%topology%nodes(p)%shape
          do z=1,d(3)
          do y=1,d(2)
          do x=1,d(1)
            if(old%nodes(p)%geometry%volume_fraction(x,y,z)<=0.or.covered(tree%topology,p,[x,y,z])) cycle
            low=old%nodes(p)%state(:,x,y,z)
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
      safe=valid
      if(.not.safe) return
      admitted_mean_rate=low_rate;admitted_mean_faces=low_flux
      if(present(redistribution_transfer)) admitted_mean_transfer=low_transfer
      ! Validate the same arithmetic that the reactor will consume. Encoding
      ! the selected mean as 2*mean-old_stage and averaging it back can erase
      ! the last trace-species inventory. Keep the mean authoritative, and
      ! retreat the entire conservative family together if rounding exhausts
      ! a budget. The already admitted low family remains the final anchor.
      do retreat=0,65
        do p=1,size(tree%nodes)
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
        do p=1,size(tree%nodes)
          if(.not.active(p).and..not.present(redistribution_transfer)) cycle
          d=tree%topology%nodes(p)%shape
          do z=1,d(3)
          do y=1,d(2)
          do x=1,d(1)
            if(old%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
            if(active(p)) then
              high=old%nodes(p)%state(:,x,y,z)+dt*(fixed_rate%fields(p)%values(:,x,y,z)+ &
                admitted_mean_rate%fields(p)%values(:,x,y,z))
            else
              if(covered(tree%topology,p,[x,y,z])) cycle
              high=old%nodes(p)%state(:,x,y,z)+admitted_mean_transfer%fields(p)%values(:,x,y,z)
            endif
            valid=valid.and.all(ieee_is_finite(high)).and.high(1)>0.and.all(high(6:)>=0)
          enddo
          enddo
          enddo
        enddo
        safe=valid
        if(safe) exit
        if(blend==0) return
        if(retreat<64) then
          blend=ieee_next_after(blend,0._dp)
        else
          blend=0
        endif
      enddo
      if(.not.safe) return
      do p=1,size(tree%nodes)
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
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        if(allocated(admitted_mean_rate%fields)) then
          rates%fields(p)%values=.5_dp*(e0%fields(p)%values+e1%fields(p)%values)
          if(present(auxiliary_source)) rates%fields(p)%values=rates%fields(p)%values+ &
            .5_dp*(a0%fields(p)%values+a1%fields(p)%values)
          rates%fields(p)%values=rates%fields(p)%values+admitted_mean_rate%fields(p)%values
        else
          rates%fields(p)%values=.5_dp*(e0%fields(p)%values+e1%fields(p)%values+ &
            d0%fields(p)%values+d1%fields(p)%values)
          if(present(auxiliary_source)) rates%fields(p)%values=rates%fields(p)%values+ &
            .5_dp*(a0%fields(p)%values+a1%fields(p)%values)
        endif
        if(.not.with_reaction.and.held_prepared.and.it==1.and.present(energy_correction)) &
          rates%fields(p)%values(5,:,:,:)=rates%fields(p)%values(5,:,:,:)+.5_dp*predictor_reset%fields(p)%values(5,:,:,:)/dt
        if(.not.with_reaction) cycle
        rates%fields(p)%values(6:,:,:,:)=rates%fields(p)%values(6:,:,:,:)+ir%fields(p)%values(:nv-5,:,:,:)
        rates%fields(p)%values(5,:,:,:)=rates%fields(p)%values(5,:,:,:)+ir%fields(p)%values(nv-4,:,:,:)
      end do
      if(.not.with_reaction.and.held_prepared.and.it==1.and.present(energy_correction)) &
        accepted_predictor_reset=.5_dp*predictor_reset_change
    end subroutine combined_rates
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
      integer :: p,x,y,z,d(3)
      real(dp) :: old_value,increment
      good=.true.;if(.not.held_auxiliary) return
      good=.false.
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        d=tree%topology%nodes(p)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          old_value=rate%fields(p)%values(1,x,y,z);increment=packet%fields(p)%values(naux,x,y,z)
          if(.not.all(ieee_is_finite([old_value,increment]))) return
          if(increment>0) then
            if(old_value>huge(1._dp)-increment) return
          else if(increment<0) then
            if(old_value< -huge(1._dp)-increment) return
          endif
          rate%fields(p)%values(1,x,y,z)=old_value+increment
        enddo
        enddo
        enddo
      enddo
      good=.true.
    end subroutine add_held_auxiliary_rate

    subroutine sources(state,external_rate,transport_rate,auxiliary_rate,stage,flux,exchange,good,fractions)
      type(amr_tree_state_3d) :: mol_transport_source,mol_transport_old
      real(dp), allocatable :: mol_body_exchange(:),mol_body_frequency
      type(amr_tree_state_3d), allocatable :: transport_held
      type(eb_reactive_tree_3d), intent(in) :: state
      type(amr_tree_state_3d), intent(out) :: external_rate,transport_rate,auxiliary_rate,exchange
      integer, intent(in) :: stage
      type(eb_integrated_flux_3d), allocatable, intent(out) :: flux(:)
      logical, intent(out) :: good
      real(dp), optional, intent(in) :: fractions(:)
      integer :: p,d(3),x,y,z
      real(dp) :: position(3),rotation_source(3)
      type(amr_tree_state_3d) :: packet,packed_state,source_gas,source_held
      type(amr_tree_state_3d) :: extra,callback_halo,callback_old
      type(eb_reactive_tree_3d) :: source_state
      type(eb_geometry_3d), allocatable :: callback_geometry(:)
      integer :: callback_padding
      type(eb_integrated_flux_3d), allocatable :: extra_flux(:)
      good=.false.
      godunov_source_stage=stage
      allocate(external_rate%fields(size(tree%nodes)))
      do p=1,size(tree%nodes)
        d=tree%topology%nodes(p)%shape
        allocate(external_rate%fields(p)%values(nv,d(1),d(2),d(3)),source=0._dp)
        if(.not.active(p)) cycle
        if(cached_native_source(stage)) then
          external_rate%fields(p)=native_source_replay%external_stages(stage)%fields(p)
          cycle
        endif
        if(present(replay_interval).and.(held_prepared.or.native_external_frozen)) then
          if(held_prepared) then
            external_rate%fields(p)=replay_interval%prepared_stages(stage)%external_rate%fields(p)
          else
            external_rate%fields(p)=replay_interval%external_stages(stage)%fields(p)
          endif
          cycle
        endif
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          if(state%nodes(p)%geometry%volume_fraction(x,y,z)<=0) cycle
          if(state%nodes(p)%state(1,x,y,z)<=0) return
          external_rate%fields(p)%values(2:4,x,y,z)=force
          external_rate%fields(p)%values(5,x,y,z)=dot_product(force,state%nodes(p)%state(2:4,x,y,z))/ &
            state%nodes(p)%state(1,x,y,z)
          if(present(frame)) then
            associate(g=>state%nodes(p)%geometry)
            position=[g%x_lower,g%y_lower,g%z_lower]+(real([x,y,z],dp)-.5_dp)*[g%dx,g%dy,g%dz]
            end associate
            call rotating_momentum_source_3d(frame,position,state%nodes(p)%state(1,x,y,z), &
              state%nodes(p)%state(2:4,x,y,z),rotation_source,good)
            if(.not.good) return
            good=.false.
            external_rate%fields(p)%values(2:4,x,y,z)=external_rate%fields(p)%values(2:4,x,y,z)+rotation_source
            ! Rotating total energy has no rotational source term.
          end if
        end do
        end do
        end do
      end do
      if(held_prepared) held_stages(stage)%external_rate=external_rate
      if(allocated(held_history_work%external_stages)) held_history_work%external_stages(stage)=external_rate
      transport_rate=external_rate
      do p=1,size(tree%nodes)
        transport_rate%fields(p)%values=0
      end do
      if(present(auxiliary_source)) then
        auxiliary_rate=transport_rate;packet=transport_rate
        if(held_auxiliary) then
          do p=1,size(tree%nodes)
            d=tree%topology%nodes(p)%shape
            deallocate(packet%fields(p)%values)
            allocate(packet%fields(p)%values(naux,d(1),d(2),d(3)),source=0._dp)
          enddo
        endif
        if(cached_native_source(stage)) then
          packet=native_source_replay%auxiliary_stages(stage)
        else if(cached_godunov_auxiliary(stage)) then
          call copy_godunov_field_3d(godunov_replay_interval%sources(stage)%auxiliary,packet,good)
          if(.not.good) return
        else if(present(replay_interval)) then
          if(held_prepared) then
            packet=replay_interval%prepared_stages(stage)%auxiliary_rate
          else
            packet=replay_interval%auxiliary_stages(stage)
          endif
        else if(.not.(mol.and.stage>1)) then
          if(held_auxiliary) then
            allocate(source_gas%fields(size(tree%nodes)))
            do p=1,size(tree%nodes)
              source_gas%fields(p)%values=state%nodes(p)%state
            enddo
            if(stage==0) then
              source_held=held_old
            else
              source_held=held_candidate
            endif
            call pack_eb_held_auxiliary_state_3d(source_gas,source_held,packed_state,good)
            if(.not.good) return
            source_state=state
            do p=1,size(tree%nodes)
              source_state%nodes(p)%state=packed_state%fields(p)%values
            enddo
            godunov_packet_pending=.true.
            call auxiliary_source(source_state,dt,stage,selected_depth,packet,good)
          else
            godunov_packet_pending=.true.
            call auxiliary_source(state,dt,stage,selected_depth,packet,good)
          endif
          if(.not.good) return
        endif
        good=.false.
        if(.not.allocated(packet%fields)) return
        if(size(packet%fields)/=size(tree%nodes).or.lbound(packet%fields,1)/=1) return
        do p=1,size(tree%nodes)
          d=tree%topology%nodes(p)%shape
          if(.not.allocated(packet%fields(p)%values)) return
          if(any(shape(packet%fields(p)%values)/=[naux,d]).or.any(lbound(packet%fields(p)%values)/=1)) return
          if(.not.active(p)) then
            packet%fields(p)%values=0;cycle
          endif
          do z=1,d(3)
          do y=1,d(2)
          do x=1,d(1)
            if(state%nodes(p)%geometry%volume_fraction(x,y,z)<=0) packet%fields(p)%values(:,x,y,z)=0
          enddo
          enddo
          enddo
          if(.not.all(ieee_is_finite(packet%fields(p)%values))) return
          auxiliary_rate%fields(p)%values=packet%fields(p)%values(:nv,:,:,:)
        enddo
        if(allocated(godunov_history)) then
          call capture_auxiliary_packet(stage,packet,good)
          if(.not.good) return
        endif
        godunov_packet_pending=.false.
        if(native_source_capture.and..not.held_prepared) &
          held_history_work%auxiliary_stages(stage)=packet
        if(stage==0) then
          auxiliary_packed0=packet
        else
          auxiliary_packed1=packet
        endif
      endif
      if(held_prepared) then
        held_stages(stage)%auxiliary_rate=transport_rate
        if(auxiliary_present) held_stages(stage)%auxiliary_rate=packet
      endif
      if(native_prefix_agreed) native_prefix_observed=max(native_prefix_observed,stage+1)
      if(prepared_redist) then
        if(prepared_godunov) then
          call prepared_diffusion_sources(state,stage,transport_rate,flux,stage_theta,good,fractions)
        else
          call prepared_mol_sources(state,stage,transport_rate,flux,stage_theta,good,fractions)
        endif
        if(.not.good) return
        if(held_auxiliary.and..not.held_prepared) then
          call add_held_auxiliary_stage(stage,good)
          if(.not.good) return
        endif
        theta=min(theta,stage_theta)
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
          if(wall_active) allocate(mol_body_exchange(nv),mol_body_frequency)
          call compute_sdc_eb_tree_transport_halos_3d(species,transport,state,extended_geometry,dt,transport_flags, &
            transport_rate,flux,mol_transport_source,mol_transport_old,stage_theta,good, &
            boundary_faces,inflow_state,inflow_temperature,wall_temperature,target_depth,ancestor_start,ancestor_end, &
            fractions,frame,embedded_wall_temperature,embedded_noslip,mol_body_exchange,mol_body_frequency, &
            inflow_mach,inflow_pressure,constant_coefficients=constant_coefficients,internal_energy=transport_held, &
            ancestor_internal_energy_start=ancestor_internal_energy_start, &
            ancestor_internal_energy_end=ancestor_internal_energy_end,coarse_linear=mol.and.present(internal_energy), &
        separate_momentum=selected_momentum, &
        consistent_velocity=selected_velocity,paired_pressure=selected_pressure,connected_directions=selected_connected, &
          primitive_from_total_energy=mol.and..not.present(frame).and..not.selected_momentum)
          if(.not.good) return
          if(wall_active) then
            heat_stage=mol_body_exchange(5);heat_frequency=mol_body_frequency
            if(noslip) wall_impulse=mol_body_exchange(2:4)
          endif
        else
        call compute_sdc_eb_tree_transport_3d(species,transport,dt,transport_flags(1),transport_flags(2), &
          transport_flags(3),transport_flags(4),state,transport_rate,flux,stage_theta,good, &
          inflow_mach,inflow_temperature,inflow_pressure,boundary_faces,inflow_state,wall_temperature, &
          target_depth,ancestor_start,ancestor_end,fractions, &
          embedded_wall_temperature=embedded_wall_temperature,extended_geometry=heat_geometry, &
          embedded_wall_energy=heat_stage,embedded_frequency=heat_frequency, &
          embedded_noslip=embedded_noslip,embedded_wall_impulse=wall_impulse,frame=frame, &
          constant_coefficients=constant_coefficients,internal_energy=transport_held, &
          ancestor_internal_energy_start=ancestor_internal_energy_start, &
          ancestor_internal_energy_end=ancestor_internal_energy_end)
        endif
        if(.not.good) return
        if(wall_active) then
          if(stage==0) then
            heat_old=heat_stage
            if(noslip) impulse_old=wall_impulse
          else
            heat_new=heat_stage
            if(noslip) impulse_new=wall_impulse
          end if
        end if
      else
        stage_theta=1
      end if
      if(mol) then
        if(present(internal_energy)) then
          if(stage==0) then
            call held_mol_source(stage,state,held_old,held_rate0,held_flux0,held_registers0,good,fractions, &
              extra,extra_flux)
          else
            call held_mol_source(stage,state,held_candidate,held_rate1,held_flux1,held_registers1,good,fractions, &
              extra,extra_flux)
          endif
        else
        call compute_mol_eb_tree_hydro_3d(species,state,mol_geometry,dt,2,extra,extra_flux,good, &
          boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
          target_depth,ancestor_start,ancestor_end,fractions,mol_lax_friedrichs,frame)
        endif
        if(.not.good) return
        do p=1,size(tree%nodes)
          if(.not.active(p)) cycle
          transport_rate%fields(p)%values=transport_rate%fields(p)%values+extra%fields(p)%values
          if(diffusive) then
            flux(p)%x=flux(p)%x+extra_flux(p)%x
            flux(p)%y=flux(p)%y+extra_flux(p)%y
            flux(p)%z=flux(p)%z+extra_flux(p)%z
          end if
        end do
        if(.not.diffusive) call move_alloc(extra_flux,flux)
        do p=1,size(tree%nodes)
          extra%fields(p)%values=0
          if(.not.active(p)) cycle
          call wall_source(state,p,extra%fields(p)%values,good,stage)
          if(.not.good) return
        end do
        if(stage==0) then
          call field_inventory(extra,mol_wall_old,good)
          mol_wall_old=dt*mol_wall_old
        else
          call field_inventory(extra,mol_wall_new,good)
          mol_wall_new=dt*mol_wall_new
        end if
        if(.not.good) return
      end if
      if(present(evaluate_spatial_diffusion).and..not.(external_diffusion_source.and.stage>1)) then
        if(present(internal_energy).and..not.allocated(transport_held)) then
          allocate(transport_held)
          if(stage==0) then
            transport_held=held_old
          else
            transport_held=held_candidate
          endif
        endif
        allocate(callback_geometry(size(tree%nodes)))
        callback_halo=transport_rate;callback_old=transport_rate
        do p=1,size(tree%nodes)
          d=tree%topology%nodes(p)%shape
          callback_padding=(extended_geometry(p)%nx-d(1))/2
          call extract_eb_geometry_3d(extended_geometry(p),[callback_padding-5,callback_padding-5,callback_padding-5], &
            d+callback_padding+6,callback_geometry(p),good)
          if(.not.good) return
          deallocate(callback_halo%fields(p)%values,callback_old%fields(p)%values)
          if(.not.active(p)) cycle
          allocate(callback_halo%fields(p)%values(nv,d(1)+12,d(2)+12,d(3)+12),source=0._dp)
        enddo
        call add_spatial_diffusion(state,callback_geometry,transport_rate,callback_halo,callback_old,flux, &
          stage_theta,good,fractions,transport_held)
        if(.not.good) return
      endif
      if(present(internal_energy).and..not.mol) then
        call initialize_held_godunov_stage(state,stage,good,fractions)
        if(.not.good) return
      endif
      if(held_auxiliary) then
        call add_held_auxiliary_stage(stage,good)
        if(.not.good) return
      endif
      if(mol.and.chemistry.and.redistribution_method=='flux_redist') then
        call capture_nonface_rate(transport_rate,flux,stage,good)
        if(.not.good) return
      endif
      if(diffusive.or.mol) then
        if(redistribution_method/='none') then
          if(redistribution_method=='state_redist'.and.redist_order>0) then
            call redist_boundary_samples(state,good)
            if(.not.good) return
          endif
          if(selected_depth>=0) then
            if(redistribution_method=='state_redist') then
              call redistribute_sdc_eb_tree_state_level_rate_3d(state,transport_rate,selected_depth,dt,redistributed, &
                exchange,good,redist_target,boundary_faces, &
        max_order=redist_order,boundary_halos=redist_halos,boundary_geometry=redist_geometry)
            else
              call redistribute_sdc_eb_tree_level_rate_3d(state,transport_rate,selected_depth,dt,redistributed, &
                exchange,good,boundary_faces)
            end if
          else
            if(redistribution_method=='state_redist') then
              call redistribute_sdc_eb_tree_state_rate_3d(state,transport_rate,dt,redistributed,good, &
                redist_target,boundary_faces, &
        max_order=redist_order,boundary_halos=redist_halos,boundary_geometry=redist_geometry)
            else
              call redistribute_sdc_eb_tree_flux_rate_3d(state,transport_rate,redistributed,good,boundary_faces)
            end if
          end if
          if(.not.good) return
          transport_rate=redistributed
        end if
        do p=1,size(tree%nodes)
          if(.not.active(p)) cycle
          call post_eb_spatial_density_rate_3d(state%nodes(p)%geometry,transport_rate%fields(p)%values,good)
          if(.not.good) return
        enddo
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
    end subroutine sources
    subroutine add_spatial_diffusion(state,padded,output,source,old_halo,flux,limiter,good,fractions,stored)
      type(amr_tree_state_3d), optional, intent(in) :: stored
      type(eb_reactive_tree_3d), intent(in) :: state
      type(eb_geometry_3d), intent(in) :: padded(:)
      type(amr_tree_state_3d), intent(inout) :: output,source,old_halo
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
      real(dp), intent(inout) :: limiter
      logical, intent(out) :: good
      real(dp), optional, intent(in) :: fractions(:)
      type(amr_tree_state_3d) :: callback_old,callback_rate,callback_source
      type(eb_integrated_flux_3d), allocatable :: callback_flux(:)
      real(dp), allocatable :: gt(:,:,:),vf(:,:,:),metrics(:,:,:),ghost_stored(:,:,:)
      real(dp) :: callback_theta
      integer(int64), allocatable :: origins(:,:),domains(:,:)
      integer, allocatable :: ratios(:)
      integer(int64) :: cell(3)
      integer :: p,q,d(3),a,x,y,z,level,index(3),width
      logical :: inside
      good=.true.
      if(.not.present(evaluate_spatial_diffusion)) return
      good=.false.;width=callback_width
      allocate(callback_old%fields(size(state%nodes)))
      callback_rate=output
      do p=1,size(state%nodes)
        d=state%topology%nodes(p)%shape
        callback_rate%fields(p)%values=0
        if(.not.active(p)) cycle
        if(any(d>huge(1)-2*width)) return
        if(width>6) then
          ! Like the MPI wide-filter route, additional donors require regular
          ! geometry. Their metrics are one; no cut-cell extrapolation occurs.
          if(any(state%nodes(p)%geometry%volume_fraction/=1).or.any(padded(p)%volume_fraction/=1)) return
          allocate(metrics(d(1)+2*width,d(2)+2*width,d(3)+2*width),source=1._dp)
        else
          metrics=padded(p)%volume_fraction
        endif
        allocate(callback_old%fields(p)%values(nv,1-width:d(1)+width,1-width:d(2)+width,1-width:d(3)+width), &
          gt(1-width:d(1)+width,1-width:d(2)+width,1-width:d(3)+width), &
          vf(1-width:d(1)+width,1-width:d(2)+width,1-width:d(3)+width))
        if(present(stored)) allocate(ghost_stored(1-width:d(1)+width,1-width:d(2)+width,1-width:d(3)+width))
        call fill_reactive_eb_tree_halo_3d(species,state,p,callback_old%fields(p)%values,gt,good, &
          boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_temperature=inflow_temperature, &
          inflow_mach=inflow_mach,inflow_pressure=inflow_pressure,ghost_volume_fraction=vf,halo_width=width, &
          target_volume_fraction=metrics,restrict_snapshot=.false.,ancestor_start=ancestor_start, &
          ancestor_end=ancestor_end,ancestor_fraction=fractions,frame=frame, &
          internal_energy=stored,ghost_internal_energy=ghost_stored, &
          ancestor_internal_energy_start=ancestor_internal_energy_start, &
          ancestor_internal_energy_end=ancestor_internal_energy_end)
        if(.not.good) return
        good=.false.
        if(allocated(old_halo%fields(p)%values)) then
          if(any(callback_old%fields(p)%values(:,-5:d(1)+6,-5:d(2)+6,-5:d(3)+6)/= &
              old_halo%fields(p)%values)) return
        else
          old_halo%fields(p)%values=callback_old%fields(p)%values(:,-5:d(1)+6,-5:d(2)+6,-5:d(3)+6)
        endif
        deallocate(gt,vf,metrics)
        if(allocated(ghost_stored)) deallocate(ghost_stored)
      enddo
      callback_theta=1
      if(cached_godunov_provider(godunov_source_stage).and.external_diffusion_source) then
        call copy_godunov_field_3d(godunov_replay_interval%sources(godunov_source_stage)%provider_rate,callback_rate,good)
        if(.not.good) return
        call copy_godunov_field_3d(godunov_replay_interval%sources(godunov_source_stage)%provider_halos,callback_source,good)
        if(.not.good) return
        call copy_godunov_fluxes_3d(godunov_replay_interval%sources(godunov_source_stage)%provider_faces,callback_flux,good)
        if(.not.good) return
        callback_theta=godunov_replay_interval%sources(godunov_source_stage)%provider_limiter
      else if(present(replay_interval).and.external_diffusion_source.and..not.held_prepared) then
        associate(saved=>replay_interval%provider_stages(1,godunov_source_stage))
        call copy_godunov_field_3d(saved%provider_rate,callback_rate,good);if(.not.good) return
        call copy_godunov_field_3d(saved%provider_halos,callback_source,good);if(.not.good) return
        call copy_godunov_fluxes_3d(saved%provider_faces,callback_flux,good);if(.not.good) return
        callback_theta=saved%provider_limiter
        end associate
      else
      godunov_packet_pending=.true.
      call evaluate_spatial_diffusion(state,padded,callback_old,dt,selected_depth, &
        callback_rate,callback_source,callback_flux,callback_theta,good,frame)
      if(.not.good) then
        prepared_provider_refused=.true.
        godunov_recovery_forbidden=.true.
        if(present(provider_refused)) provider_refused=.true.
        return
      endif
      endif
      good=.false.
      if(.not.allocated(callback_rate%fields).or..not.allocated(callback_source%fields).or. &
        .not.allocated(callback_flux)) return
      if(size(callback_rate%fields)/=size(state%nodes).or.size(callback_source%fields)/=size(state%nodes).or. &
        size(callback_flux)/=size(state%nodes)) return
      if(lbound(callback_rate%fields,1)/=1.or.lbound(callback_source%fields,1)/=1.or.lbound(callback_flux,1)/=1) return
      if(.not.ieee_is_finite(callback_theta).or.callback_theta<0.or.callback_theta>1) return
      do p=1,size(state%nodes)
        d=state%topology%nodes(p)%shape
        if(.not.allocated(callback_rate%fields(p)%values)) return
        if(any(shape(callback_rate%fields(p)%values)/=[nv,d]).or. &
          any(lbound(callback_rate%fields(p)%values)/=1)) return
        if(.not.all(ieee_is_finite(callback_rate%fields(p)%values))) return
        if(.not.active(p)) then
          if(allocated(callback_source%fields(p)%values)) return
          if(any(callback_rate%fields(p)%values/=0)) return
          cycle
        endif
        if(.not.allocated(callback_source%fields(p)%values)) return
        if(any(shape(callback_source%fields(p)%values)/=[nv,d+12]).or. &
          any(lbound(callback_source%fields(p)%values)/=1)) return
        if(.not.all(ieee_is_finite(callback_source%fields(p)%values))) return
        if(any(callback_source%fields(p)%values(:,7:d(1)+6,7:d(2)+6,7:d(3)+6)/= &
          callback_rate%fields(p)%values)) return
        if(.not.valid_eb_integrated_flux_3d(callback_flux(p),state%nodes(p)%geometry,nv)) return
      enddo
      if(allocated(godunov_history)) then
        call capture_provider_packet(godunov_source_stage,1,callback_rate,callback_source,callback_flux,callback_theta,good)
        if(.not.good) return
      endif
      if(allocated(held_history_work%provider_stages)) then
        associate(saved=>held_history_work%provider_stages(1,godunov_source_stage))
        call copy_godunov_field_3d(callback_rate,saved%provider_rate,good);if(.not.good) return
        call copy_godunov_field_3d(callback_source,saved%provider_halos,good);if(.not.good) return
        call copy_godunov_fluxes_3d(callback_flux,saved%provider_faces,good);if(.not.good) return
        saved%provider_accepted=.true.;saved%provider_kind=1;saved%index=godunov_source_stage
        saved%provider_limiter=callback_theta
        end associate
      endif
      godunov_packet_pending=.false.
      do p=1,size(state%nodes)
        if(.not.active(p)) cycle
        output%fields(p)%values=output%fields(p)%values+callback_rate%fields(p)%values
        source%fields(p)%values=source%fields(p)%values+callback_source%fields(p)%values
        flux(p)%x=flux(p)%x+callback_flux(p)%x
        flux(p)%y=flux(p)%y+callback_flux(p)%y
        flux(p)%z=flux(p)%z+callback_flux(p)%z
        if(.not.all(ieee_is_finite(output%fields(p)%values)).or. &
          .not.all(ieee_is_finite(source%fields(p)%values))) return
        if(.not.valid_eb_integrated_flux_3d(flux(p),state%nodes(p)%geometry,nv)) return
      enddo
      ! Recipient valid cells and periodic images retain the actual combined
      ! stage rate before StateRedist evaluates neighborhoods/registers.
      call amr_tree_integer_layout_3d(state%topology,origins,domains,ratios,good)
      if(.not.good) return
      good=.false.
      do p=1,size(state%nodes)
        if(.not.active(p)) cycle
        d=state%topology%nodes(p)%shape;level=state%topology%nodes(p)%level
        do z=4,d(3)+9
        do y=4,d(2)+9
        do x=4,d(1)+9
          cell=origins(:,p)+int([x,y,z]-7,int64);inside=.true.
          do a=1,3
            if(boundary_faces(2*a-1)=='periodic') then
              cell(a)=modulo(cell(a),domains(a,level+1))
            else if(cell(a)<0.or.cell(a)>=domains(a,level+1)) then
              inside=.false.
            endif
          enddo
          if(.not.inside) cycle
          do q=1,size(state%nodes)
            if(state%topology%nodes(q)%level/=level) cycle
            if(any(cell<origins(:,q)).or.any(cell>=origins(:,q)+int(state%topology%nodes(q)%shape,int64))) cycle
            index=int(cell-origins(:,q))+1
            source%fields(p)%values(:,x,y,z)=output%fields(q)%values(:,index(1),index(2),index(3))
            exit
          enddo
        enddo
        enddo
        enddo
      enddo
      limiter=min(limiter,callback_theta);good=.true.
    end subroutine add_spatial_diffusion
    subroutine prepared_diffusion_sources(state,stage,output,flux,limiter,good,fractions)
      type(eb_reactive_tree_3d), intent(in) :: state
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
      integer :: p,m(3),depth,first_depth,last_depth
      if(present(internal_energy)) then
        allocate(transport_held)
        if(stage==0) then
          transport_held=held_old
        else
          transport_held=held_candidate
        endif
      endif
      if(wall_active) allocate(body_exchange(nv),frequency)
      if(diffusive) then
        if(present(evaluate_spatial_diffusion)) then
        call compute_sdc_eb_tree_transport_halos_3d(species,transport,state,register_geometry,dt,transport_flags, &
          output,flux,source,old_halo,limiter,good,boundary_faces,inflow_state,inflow_temperature,wall_temperature, &
          target_depth,ancestor_start,ancestor_end,fractions,frame,embedded_wall_temperature,embedded_noslip, &
          body_exchange,frequency,inflow_mach,inflow_pressure, &
          constant_coefficients=constant_coefficients,internal_energy=transport_held, &
          ancestor_internal_energy_start=ancestor_internal_energy_start, &
          ancestor_internal_energy_end=ancestor_internal_energy_end)
        else
        call compute_sdc_eb_tree_transport_halos_3d(species,transport,state,register_geometry,dt,transport_flags, &
          output,flux,source,old_halo,limiter,good,boundary_faces,inflow_state,inflow_temperature,wall_temperature, &
          target_depth,ancestor_start,ancestor_end,fractions,frame,embedded_wall_temperature,embedded_noslip, &
          body_exchange,frequency,inflow_mach,inflow_pressure,registers,redist_order,redist_target, &
          constant_coefficients=constant_coefficients,internal_energy=transport_held, &
          ancestor_internal_energy_start=ancestor_internal_energy_start, &
          ancestor_internal_energy_end=ancestor_internal_energy_end)
        endif
      else
        call compute_sdc_eb_tree_transport_halos_3d(species,transport,state,register_geometry,dt,transport_flags, &
          output,flux,source,old_halo,limiter,good,boundary_faces,inflow_state,inflow_temperature,wall_temperature, &
          target_depth,ancestor_start,ancestor_end,fractions, &
          constant_coefficients=constant_coefficients,internal_energy=transport_held, &
          ancestor_internal_energy_start=ancestor_internal_energy_start, &
          ancestor_internal_energy_end=ancestor_internal_energy_end)
        if(.not.good) return
        allocate(registers(size(tree%nodes)))
        do p=1,size(tree%nodes)
          if(.not.active(p)) cycle
          call initialize_state_redist_tree_register_3d(tree%topology,p,nv,6,registers(p),good,boundary_faces,1._dp)
          if(.not.good) return
        enddo
      endif
      if(.not.good) return
      if(present(evaluate_spatial_diffusion)) then
        call add_spatial_diffusion(state,register_geometry,output,source,old_halo,flux,limiter,good,fractions,transport_held)
        if(.not.good) return
        first_depth=max(0,selected_depth);last_depth=first_depth
        if(selected_depth<0) last_depth=tree%topology%level_count()-1
        do depth=first_depth,last_depth
          call compute_state_redist_level_stage_3d(state%topology,register_geometry,old_halo,source,nv, &
            depth,dt,output,registers,good,boundary_faces,redist_target,redist_order,register_scale=1._dp)
          if(.not.good) return
        enddo
      endif
      if(wall_active) then
        if(stage==0) then
          heat_old=body_exchange(5);impulse_old=body_exchange(2:4)
        else
          heat_new=body_exchange(5);impulse_new=body_exchange(2:4)
        endif
      endif
      ! Post-redistribution density consistency precedes gas/EOS stepping.
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        call post_eb_spatial_density_rate_3d(state%nodes(p)%geometry,output%fields(p)%values,good)
        if(.not.good) return
      enddo
      if(present(internal_energy)) then
        call initialize_held_godunov_stage(state,stage,good,fractions)
        if(.not.good) return
      endif
      if(stage==0) then
        call move_alloc(registers,stage_registers0)
      else
        if(stage>1) then
          do p=1,size(tree%nodes)
            if(.not.active(p)) cycle
            m=tree%topology%nodes(p)%shape
            call accumulate_state_redist_coarse_register_3d(registers(p),stage_registers1(p),nv,m+12,good)
            if(.not.good) return
          enddo
        endif
        call move_alloc(registers,stage_registers1)
      endif
    end subroutine prepared_diffusion_sources
    subroutine prepared_mol_sources(state,stage,output,flux,limiter,good,fractions)
      type(amr_tree_state_3d), allocatable :: transport_held
      type(eb_reactive_tree_3d), intent(in) :: state
      integer, intent(in) :: stage
      type(amr_tree_state_3d), intent(inout) :: output
      type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
      real(dp), intent(inout) :: limiter
      logical, intent(out) :: good
      real(dp), optional, intent(in) :: fractions(:)
      type(amr_tree_state_3d) :: source,old_halo,mol_source,mol_old,extra,scalar_raw,scalar_old
      type(eb_integrated_flux_3d), allocatable :: mol_flux(:)
      type(state_redist_register_3d), allocatable :: registers(:)
      real(dp), allocatable :: body_exchange(:),frequency
      integer :: p,depth,first_depth,last_depth,m(3)
      good=.false.
      if(wall_active) allocate(body_exchange(nv),frequency)
      if(present(internal_energy)) then
        allocate(transport_held)
        if(stage==0) then
          transport_held=held_old
        else
          transport_held=held_candidate
        endif
      endif
      call compute_sdc_eb_tree_transport_halos_3d(species,transport,state,extended_geometry,dt,transport_flags, &
        output,flux,source,old_halo,limiter,good,boundary_faces,inflow_state,inflow_temperature,wall_temperature, &
        target_depth,ancestor_start,ancestor_end,fractions,frame,embedded_wall_temperature,embedded_noslip, &
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
        if(stage==0) then
          call held_mol_source(stage,state,held_old,held_rate0,held_flux0,held_registers0,good,fractions, &
            extra,mol_flux,mol_source,mol_old,scalar_raw,scalar_old)
        else
          call held_mol_source(stage,state,held_candidate,held_rate1,held_flux1,held_registers1,good,fractions, &
            extra,mol_flux,mol_source,mol_old,scalar_raw,scalar_old)
        endif
      else
      call compute_mol_eb_tree_hydro_3d(species,state,extended_geometry,dt,2,extra,mol_flux,good, &
        boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
        target_depth,ancestor_start,ancestor_end,fractions,mol_lax_friedrichs,frame, &
        spatial_rate_halos=mol_source,old_state_halos=mol_old)
      endif
      if(.not.good) return
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        ! Both operators evaluated the same prepared old state at this level.
        if(any(old_halo%fields(p)%values/=mol_old%fields(p)%values)) then
          good=.false.;return
        endif
        output%fields(p)%values=output%fields(p)%values+extra%fields(p)%values
        source%fields(p)%values=source%fields(p)%values+mol_source%fields(p)%values
        flux(p)%x=flux(p)%x+mol_flux(p)%x
        flux(p)%y=flux(p)%y+mol_flux(p)%y
        flux(p)%z=flux(p)%z+mol_flux(p)%z
      enddo
      if(.not.(external_diffusion_source.and.stage>1)) then
        call add_spatial_diffusion(state,extended_geometry,output,source,old_halo,flux,limiter,good,fractions, &
          transport_held)
        if(.not.good) return
      endif
      ! The pinned Diffusion MOL factor is a ledger weight, never a RHS weight.
      first_depth=max(0,selected_depth);last_depth=first_depth
      if(selected_depth<0) last_depth=tree%topology%level_count()-1
      do depth=first_depth,last_depth
        call compute_state_redist_level_stage_3d(state%topology,extended_geometry,old_halo,source,nv, &
          depth,dt,output,registers,good,boundary_faces,redist_target,redist_order,register_scale=.5_dp)
        if(.not.good) return
      enddo
      ! Post-redistribution density consistency precedes gas/EOS stepping.
      do p=1,size(tree%nodes)
        if(.not.active(p)) cycle
        call post_eb_spatial_density_rate_3d(state%nodes(p)%geometry,output%fields(p)%values,good)
        if(.not.good) return
      enddo
      if(held_prepared) then
        held_stages(stage)%source_stage=stage;held_stages(stage)%depth=selected_depth
        held_stages(stage)%interval=dt;held_stages(stage)%max_order=redist_order
        held_stages(stage)%target_volume_fraction=redist_target
        held_stages(stage)%state_halos=old_halo;held_stages(stage)%scalar_state_halos=scalar_old
        held_stages(stage)%raw=source;held_stages(stage)%rate=output
        held_stages(stage)%faces=flux;held_stages(stage)%registers=registers
        held_stages(stage)%scalar_raw=scalar_raw
        if(stage==0) then
          held_stages(stage)%scalar_rate=held_rate0;held_stages(stage)%scalar_faces=held_flux0
          held_stages(stage)%scalar_registers=held_registers0
          call blend_held_prepared_stage(stage,held_spatial_weight,output,flux,registers, &
            held_rate0,held_flux0,held_registers0,good)
        else
          held_stages(stage)%scalar_rate=held_rate1;held_stages(stage)%scalar_faces=held_flux1
          held_stages(stage)%scalar_registers=held_registers1
          call blend_held_prepared_stage(stage,held_spatial_weight,output,flux,registers, &
            held_rate1,held_flux1,held_registers1,good)
        endif
        if(.not.good) return
      endif
      if(stage==0) then
        call move_alloc(registers,stage_registers0)
      else
        ! Coarse storage is additive on every fixed-reference evaluation,
        ! including provisional endpoints with zero accepted flux weight.
        ! Fine dm follows the latest accepted endpoint flux only.
        if(stage>1.and..not.(chemistry.and.present(internal_energy))) then
          do p=1,size(tree%nodes)
            if(.not.active(p)) cycle
            m=tree%topology%nodes(p)%shape
            call accumulate_state_redist_coarse_register_3d(registers(p),stage_registers1(p),nv,m+12,good)
            if(.not.good) return
          enddo
        endif
        call move_alloc(registers,stage_registers1)
      endif
      do p=1,size(tree%nodes)
        extra%fields(p)%values=0
        if(.not.active(p)) cycle
        call wall_source(state,p,extra%fields(p)%values,good,stage)
        if(.not.good) return
      enddo
      if(stage==0) then
        call field_inventory(extra,mol_wall_old,good)
        mol_wall_old=dt*mol_wall_old
      else
        call field_inventory(extra,mol_wall_new,good)
        mol_wall_new=dt*mol_wall_new
      endif
    end subroutine prepared_mol_sources
    subroutine wall_source(state,p,output,good,stage)
      type(eb_reactive_tree_3d), intent(in) :: state
      integer, intent(in) :: p,stage
      real(dp), intent(out) :: output(:,:,:,:)
      logical, intent(out) :: good
      real(dp), allocatable :: potential(:,:,:),stored(:,:,:)
      real(dp) :: position(3)
      integer :: x,y,z
      good=.false.;output=0
      associate(g=>state%nodes(p)%geometry)
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
      call reactive_eb_slip_wall_source_3d(species,state%nodes(p)%state,state%nodes(p)%temperature, &
        g,output,good,potential,stored)
      end associate
    end subroutine wall_source
    ! Reset the selected held scalar before EOS recovery. Only exposed fluid
    ! cells contribute to the numerical-energy ledger, matching other common
    ! interval inventories; covered coarse states still receive their reset.

    subroutine recover(good,reset_phase)
      logical, intent(out) :: good
      integer, optional, intent(in) :: reset_phase
      good=.false.
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
        call average_down_reactive_eb_tree_3d(species,candidate,good,frame,held_candidate)
        if(.not.good) return
      endif
      if(present(reset_phase).and.present(energy_correction)) then
        call reset_held_eb_tree_family_3d(species,candidate,held_candidate,reset_delta,held_policy,selected_depth, &
          reset_change,good,frame, &
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
          if(.not.ieee_is_finite(closure_inventory(1))) then
            good=.false.;return
          endif
        end select
        return
      endif
      do node=1,size(tree%nodes)
        if(.not.active(node)) cycle
        n=tree%topology%nodes(node)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(candidate%nodes(node)%geometry%volume_fraction(i,j,k)<=0) cycle
          ! Native construct_Snew starts each nonreacting CTU update from
          ! S_old, including UTEMP. The pinned inverse can retain its seed.
          ! Limit this ownership correction to the total-energy pinned path.
          if(.not.mol.and..not.chemistry.and..not.present(internal_energy).and. &
             species(1)%pinned_eos_extrapolation.and.present(reset_phase)) then
            if(reset_phase==1.or.reset_phase==2) candidate%nodes(node)%temperature(i,j,k)=old%nodes(node)%temperature(i,j,k)
          endif
          call recover_reactive_eb_cell_3d(species,candidate%nodes(node)%geometry,[i,j,k],candidate%nodes(node)%state(:,i,j,k), &
            candidate%nodes(node)%temperature(i,j,k),primitive,recovered,speed,valid,frame)
          if(.not.valid) return
          candidate%nodes(node)%temperature(i,j,k)=recovered
        end do
        end do
        end do
      end do
      good=.true.
    end subroutine recover
    subroutine field_inventory(field,value,good)
      type(amr_tree_state_3d), intent(in) :: field
      real(dp), intent(inout) :: value(:)
      logical, intent(out) :: good
      type(eb_reactive_tree_3d) :: integral_tree
      integer :: p
      integral_tree=old
      do p=1,size(tree%nodes)
        integral_tree%nodes(p)%state=field%fields(p)%values
      end do
      call integrate_reactive_eb_tree_3d(integral_tree,value,good)
    end subroutine field_inventory
  end subroutine advance_reactive_eb_tree_sdc_3d

  subroutine advance_reactive_eb_tree_full_3d(species,reactions,transport,solver,dt,chemistry_enabled, &
      rtol,atol,transport_enabled,viscosity,conduction,diffusion,barodiffusion,tree,boundary_change, &
      wall_change,chemical_change,minimum_theta,ok,integrator,subcycling, &
      inflow_mach,inflow_temperature,inflow_pressure,difmag,extended_geometry,boundary_faces,inflow_state,wall_temperature)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt,rtol,atol
    logical, intent(in) :: chemistry_enabled,transport_enabled,viscosity,conduction,diffusion,barodiffusion
    type(eb_reactive_tree_3d), intent(inout) :: tree
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemical_change(:),minimum_theta
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: wall_temperature(6)
    real(dp), optional, intent(in) :: difmag
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: integrator
    logical, optional, intent(in) :: subcycling
    type(eb_reactive_tree_3d) :: candidate
    real(dp), allocatable :: boundary(:),wall(:),chemical(:),exchange(:),transport_sum(:)
    real(dp) :: theta,stage_theta
    integer :: nv,half
    logical :: valid,recursive_steps
    ok=.false.;nv=reactive_nvar(size(species))
    if (.not.ieee_is_finite(dt).or.dt<=0._dp) return
    if (size(boundary_change)/=nv.or.size(wall_change)/=nv.or.size(chemical_change)/=nv) return
    if (.not.valid_tree(tree,nv)) return
    recursive_steps=.false.
    if (present(subcycling)) recursive_steps=subcycling
    if(present(difmag)) then
      if(.not.ieee_is_finite(difmag).or.difmag<0) return
    end if
    if(present(wall_temperature)) then
      if(.not.present(boundary_faces)) return
      if(.not.all(ieee_is_finite(wall_temperature)).or.any(wall_temperature<=0)) return
    end if
    if(present(boundary_faces)) then
      if(any(boundary_faces=='isothermal_wall')) then
        if(.not.present(wall_temperature).or..not.transport_enabled.or..not.conduction) return
      end if
    end if
    candidate=tree
    allocate(boundary(nv),wall(nv),chemical(nv),exchange(nv),transport_sum(nv))
    chemical=0._dp;transport_sum=0._dp;theta=1._dp
    do half=1,2
      if (half==1.and.chemistry_enabled) then
        call reaction_half(valid)
        if (.not.valid) return
      end if
      if (transport_enabled) then
        if (recursive_steps) then
          call advance_reactive_eb_tree_subcycle_transport_3d(species,transport,0.5_dp*dt,viscosity,conduction, &
            diffusion,barodiffusion,candidate,exchange,stage_theta,valid, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
            boundary_faces=boundary_faces,inflow_state=inflow_state,wall_temperature=wall_temperature)
        else
          call advance_reactive_eb_tree_transport_3d(species,transport,0.5_dp*dt,viscosity,conduction, &
            diffusion,barodiffusion,candidate,exchange,stage_theta,valid, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
            boundary_faces=boundary_faces,inflow_state=inflow_state,wall_temperature=wall_temperature)
        end if
        if (.not.valid) return
        transport_sum=transport_sum+exchange;theta=min(theta,stage_theta)
      end if
      if (half==1) then
        if (recursive_steps) then
          call advance_reactive_eb_tree_subcycle_hydro_3d(species,solver,dt,candidate,boundary,wall,valid, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
            boundary_faces=boundary_faces,inflow_state=inflow_state, &
            difmag=difmag,extended_geometry=extended_geometry)
        else
          call advance_reactive_eb_tree_hydro_3d(species,solver,dt,candidate,boundary,wall,valid, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
            boundary_faces=boundary_faces,inflow_state=inflow_state, &
            difmag=difmag,extended_geometry=extended_geometry)
        end if
        if (.not.valid) return
      else if (chemistry_enabled) then
        call reaction_half(valid)
        if (.not.valid) return
      end if
    end do
    boundary=boundary+transport_sum
    if (.not.all(ieee_is_finite(boundary)).or..not.all(ieee_is_finite(chemical))) return
    tree=candidate;boundary_change=boundary;wall_change=wall
    chemical_change=chemical;minimum_theta=theta;ok=.true.
  contains
    subroutine reaction_half(good)
      logical, intent(out) :: good
      call react_reactive_eb_tree_3d(species,reactions,0.5_dp*dt,rtol,atol,candidate,exchange,good,integrator)
      if (good) chemical=chemical+exchange
    end subroutine reaction_half
  end subroutine advance_reactive_eb_tree_full_3d

  ! Two-cell donor margins for the quadratic embedded-wall stencil. Publish
  ! composite wall exchange only after every active patch succeeds.
  ! Transport has zero total-mass flux, so its energy flux needs no potential
  ! correction. Restore ordinary energy only in private input copies, after
  ! conserved rotating halos have been filled at their receiving positions.
  subroutine ordinary_transport_energy(frame,geometry,first,active,state,ok)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(eb_geometry_3d), intent(in) :: geometry
    integer, intent(in) :: first(3)
    logical, intent(in) :: active(:,:,:)
    real(dp), intent(inout) :: state(:,:,:,:)
    logical, intent(out) :: ok
    real(dp) :: position(3),energy
    integer :: i,j,k
    ok=.true.
    if(.not.present(frame)) return
    ok=.false.
    if(.not.frame%is_valid().or.size(state,1)<5) return
    if(any(shape(active)/=shape(state(1,:,:,:)))) return
    do k=1,size(active,3)
    do j=1,size(active,2)
    do i=1,size(active,1)
      if(.not.active(i,j,k)) cycle
      position=[geometry%x_lower,geometry%y_lower,geometry%z_lower]+ &
        (real(first+[i,j,k]-1,dp)-.5_dp)*[geometry%dx,geometry%dy,geometry%dz]
      energy=state(5,i,j,k)
      call convert_rotating_energy_3d(frame,position,state(1,i,j,k),state(5,i,j,k),.false.,energy,ok)
      if(.not.ok) return
      state(5,i,j,k)=energy
    end do
    end do
    end do
    ok=.true.
  end subroutine ordinary_transport_energy

  subroutine compute_eb_tree_wall_transport_3d(species,transport,tree,extended_geometry, &
      conduction,noslip,wall_temperature,rate,wall_exchange,frequency,ok, &
      boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,target_depth, &
      ancestor_start,ancestor_end,ancestor_fraction,frame, &
      constant_coefficients,internal_energy,ancestor_internal_energy_start,ancestor_internal_energy_end)
    real(dp), optional, intent(in) :: constant_coefficients(:)
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy,ancestor_internal_energy_start, &
      ancestor_internal_energy_end
    real(dp), allocatable :: ghost_held(:,:,:),patch_held(:,:,:)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: extended_geometry(:)
    logical, intent(in) :: conduction,noslip
    real(dp), intent(in) :: wall_temperature
    type(amr_tree_state_3d), intent(inout) :: rate
    real(dp), intent(inout) :: wall_exchange(:),frequency
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    integer, optional, intent(in) :: target_depth
    type(eb_reactive_tree_3d), optional, intent(in) :: ancestor_start,ancestor_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    type(amr_tree_state_3d) :: candidate
    type(eb_reactive_tree_3d) :: ledger
    type(eb_geometry_3d) :: center
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),vf(:,:,:)
    real(dp) :: bound,node_frequency,inventory(reactive_nvar(size(species))),node_exchange(reactive_nvar(size(species)))
    integer :: nv,node,n(3),depth
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species));depth=-1
    if(present(internal_energy)) then
      if(.not.internal_energy%is_valid(tree%topology)) return
      if(size(internal_energy%fields(1)%values,1)/=1) return
    endif
    if(present(constant_coefficients)) then
      if(.not.valid_constant_transport_coefficients(constant_coefficients)) return
    endif
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    if(.not.valid_tree(tree,nv).or.size(transport)/=size(species)) return
    if(.not.ieee_is_finite(wall_temperature).or.size(wall_exchange)/=nv) return
    if(conduction.and.wall_temperature<=0) return
    if(size(extended_geometry)/=size(tree%nodes)) return
    if(present(target_depth)) then
      depth=target_depth
      if(depth<0.or.depth>=tree%topology%level_count()) return
    end if
    if(present(ancestor_start).or.present(ancestor_end).or.present(ancestor_fraction)) then
      if(depth<0) return
    end if
    allocate(candidate%fields(size(tree%nodes)));bound=0
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      if(any([extended_geometry(node)%nx,extended_geometry(node)%ny,extended_geometry(node)%nz]/=n+4)) return
      call extract_eb_geometry_3d(extended_geometry(node),[3,3,3],n+2,center,valid)
      if(.not.valid) return
      if(.not.matching_eb_metrics_3d(center,tree%nodes(node)%geometry)) return
      allocate(candidate%fields(node)%values(nv,n(1),n(2),n(3)),source=0._dp)
      if(depth>=0.and.tree%topology%nodes(node)%level/=depth) cycle
      allocate(ghost(nv,-1:n(1)+2,-1:n(2)+2,-1:n(3)+2), &
        gt(-1:n(1)+2,-1:n(2)+2,-1:n(3)+2),vf(-1:n(1)+2,-1:n(2)+2,-1:n(3)+2))
      if(present(internal_energy)) then
        allocate(ghost_held(-1:n(1)+2,-1:n(2)+2,-1:n(3)+2))
        patch_held=internal_energy%fields(node)%values(1,:,:,:)
      endif
      call fill_reactive_eb_tree_halo_3d(species,tree,node,ghost,gt,valid, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
        inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,ghost_volume_fraction=vf, &
        halo_width=2,target_volume_fraction=extended_geometry(node)%volume_fraction,restrict_snapshot=.false., &
        ancestor_start=ancestor_start,ancestor_end=ancestor_end,ancestor_fraction=ancestor_fraction,frame=frame, &
        internal_energy=internal_energy,ghost_internal_energy=ghost_held, &
        ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end)
      if(.not.valid) return
      if(.not.present(internal_energy)) then
      call ordinary_transport_energy(frame,tree%nodes(node)%geometry,[-1,-1,-1],vf>0,ghost,valid)
      endif
      if(.not.valid) return
      call compute_eb_wall_transport_source_3d(species,transport,ghost,gt,extended_geometry(node), &
        [3,3,3],n+2,conduction,noslip,wall_temperature,candidate%fields(node)%values,node_exchange,node_frequency,valid, &
        constant_coefficients=constant_coefficients,internal_energy_density=ghost_held)
      if(.not.valid) return
      bound=max(bound,node_frequency)
      if(allocated(ghost_held)) deallocate(ghost_held)
      deallocate(ghost,gt,vf)
    end do
    ledger=tree
    do node=1,size(tree%nodes)
      ledger%nodes(node)%state=candidate%fields(node)%values
    end do
    ! Exclude covered coarse cells rather than summing per-patch ledgers.
    call integrate_reactive_eb_tree_3d(ledger,inventory,valid)
    if(.not.valid) return
    call move_alloc(candidate%fields,rate%fields)
    wall_exchange=inventory;frequency=bound;ok=.true.
  end subroutine compute_eb_tree_wall_transport_3d

  ! Common-interval molecular transport: closed outflow edges or open Mach inlets.
  ! Provisional transport rates and time-integrated face densities on every
  ! level. Parent states are retained until the coupled driver synchronizes.
  ! Raw transport source halos at the selected level, prior to StateRedist.
  ! Published old states retain rotating energy; all transport/EB wall thermo
  ! uses a private ordinary-energy copy at the actual ghost coordinates.
  subroutine compute_sdc_eb_tree_transport_halos_3d(species,transport,tree,extended_geometry,dt,flags, &
      rate,flux,source_halos,old_halos,minimum_theta,ok,boundary_faces,inflow_state,inflow_temperature, &
      wall_temperature,target_depth,ancestor_start,ancestor_end,ancestor_fraction,frame, &
      embedded_wall_temperature,embedded_noslip,embedded_exchange,embedded_frequency,inflow_mach,inflow_pressure, &
      registers,max_order,target_volume_fraction, &
      constant_coefficients,internal_energy,ancestor_internal_energy_start,ancestor_internal_energy_end,coarse_linear, &
      separate_momentum,consistent_velocity,paired_pressure,connected_directions,primitive_from_total_energy)
    logical, optional, intent(in) :: separate_momentum,consistent_velocity,paired_pressure,connected_directions
    logical, optional, intent(in) :: coarse_linear
    ! Normal face transport may share the ordinary MOL total-energy primitive recovery.
    ! Held halos and wall-source recovery retain their existing contracts.
    logical, optional, intent(in) :: primitive_from_total_energy
    logical :: linear_halos,use_total_energy
    real(dp), optional, intent(in) :: constant_coefficients(:)
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy,ancestor_internal_energy_start, &
      ancestor_internal_energy_end
    real(dp), allocatable :: ghost_held(:,:,:),patch_held(:,:,:)
    type(state_redist_register_3d), allocatable, optional, intent(inout) :: registers(:)
    integer, optional, intent(in) :: max_order
    real(dp), optional, intent(in) :: target_volume_fraction
    type(state_redist_register_3d), allocatable :: pending_registers(:)
    integer :: first_depth,last_depth,stage_depth
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: extended_geometry(:)
    real(dp), intent(in) :: dt
    logical, intent(in) :: flags(4)
    type(amr_tree_state_3d), intent(inout) :: rate,source_halos,old_halos
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    real(dp), intent(inout) :: minimum_theta
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(6),wall_temperature(6),inflow_mach(6), &
      inflow_pressure(6)
    integer, optional, intent(in) :: target_depth
    type(eb_reactive_tree_3d), optional, intent(in) :: ancestor_start,ancestor_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: embedded_wall_temperature
    logical, optional, intent(in) :: embedded_noslip
    real(dp), optional, intent(inout) :: embedded_exchange(:),embedded_frequency
    type(amr_tree_state_3d) :: rates,halos,old,body
    type(eb_reactive_tree_3d) :: ledger
    type(eb_integrated_flux_3d), allocatable :: pending_flux(:)
    type(eb_geometry_3d) :: center
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: ghost(:,:,:,:),temp(:,:,:),vf(:,:,:),ordinary(:,:,:,:),patch_halo(:,:,:,:), &
      body_halo(:,:,:,:),valid_body(:,:,:,:)
    real(dp) :: theta,patch_theta,frequency,patch_frequency,ignored_frequency,wall_t, &
      exchange(reactive_nvar(size(species))),ignored_exchange(reactive_nvar(size(species)))
    integer :: p,nv,n(3),depth,level
    logical :: good,faces(6),noslip,wall_active
    character(len=32) :: labels(6)
    ok=.false.;nv=reactive_nvar(size(species));depth=-1;theta=1;frequency=0;wall_t=0;noslip=.false.
    use_total_energy=.false.
    if(present(primitive_from_total_energy)) use_total_energy=primitive_from_total_energy
    linear_halos=.false.
    if(present(coarse_linear)) linear_halos=coarse_linear
    if(use_total_energy.and.present(frame)) return
    if(present(separate_momentum)) then
      if(use_total_energy.and.separate_momentum) return
    endif
    if(linear_halos.and..not.present(internal_energy)) return
    if(present(consistent_velocity)) then
      if(consistent_velocity) then
        if(.not.present(separate_momentum)) return
        if(.not.separate_momentum) return
      endif
    endif
    if(present(connected_directions)) then
      if(connected_directions) then
        if(.not.present(separate_momentum)) return
        if(.not.separate_momentum) return
      endif
    endif
    if(present(separate_momentum)) then
      if(separate_momentum.and..not.linear_halos) return
    endif
    if(present(internal_energy)) then
      if(.not.internal_energy%is_valid(tree%topology)) return
      if(size(internal_energy%fields(1)%values,1)/=1) return
    endif
    if(present(constant_coefficients)) then
      if(.not.valid_constant_transport_coefficients(constant_coefficients)) return
    endif
    if((present(max_order).or.present(target_volume_fraction)).and..not.present(registers)) return
    if(.not.valid_tree(tree,nv).or.size(transport)/=size(species)) return
    if(.not.ieee_is_finite(dt).or.dt<=0) return
    if(size(extended_geometry)/=size(tree%nodes)) return
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    endif
    if(present(wall_temperature)) then
      if(.not.present(boundary_faces)) return
      if(.not.all(ieee_is_finite(wall_temperature)).or.any(wall_temperature<=0)) return
    endif
    if(present(target_depth)) depth=target_depth
    if(depth< -1.or.depth>=tree%topology%level_count()) return
    if(present(target_depth).and.depth<0) return
    if(present(ancestor_start).or.present(ancestor_end).or.present(ancestor_fraction)) then
      if(depth<0) return
    endif
    if(present(embedded_noslip)) noslip=embedded_noslip
    if(noslip.and..not.flags(1)) return
    wall_active=noslip.or.present(embedded_wall_temperature)
    if(wall_active.neqv.present(embedded_exchange)) return
    if(wall_active.neqv.present(embedded_frequency)) return
    if(present(embedded_exchange)) then
      if(size(embedded_exchange)/=nv) return
    endif
    if(present(embedded_wall_temperature)) then
      wall_t=embedded_wall_temperature
      if(.not.flags(2).or..not.ieee_is_finite(wall_t).or.wall_t<=0) return
    endif
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,good)
    if(.not.good) return
    allocate(rates%fields(size(tree%nodes)),halos%fields(size(tree%nodes)),old%fields(size(tree%nodes)), &
      body%fields(size(tree%nodes)),pending_flux(size(tree%nodes)))
    do p=1,size(tree%nodes)
      n=tree%topology%nodes(p)%shape
      if(any([extended_geometry(p)%nx,extended_geometry(p)%ny,extended_geometry(p)%nz]/=n+12)) return
      call extract_eb_geometry_3d(extended_geometry(p),[7,7,7],n+6,center,good)
      if(.not.good) return
      if(.not.matching_eb_metrics_3d(center,tree%nodes(p)%geometry)) return
      allocate(rates%fields(p)%values(nv,n(1),n(2),n(3)),source=0._dp)
      allocate(body%fields(p)%values(nv,n(1),n(2),n(3)),source=0._dp)
      if(depth>=0.and.tree%topology%nodes(p)%level/=depth) cycle
      allocate(ghost(nv,1-6:n(1)+6,1-6:n(2)+6,1-6:n(3)+6), &
        temp(1-6:n(1)+6,1-6:n(2)+6,1-6:n(3)+6),vf(1-6:n(1)+6,1-6:n(2)+6,1-6:n(3)+6))
      if(present(internal_energy)) then
        allocate(ghost_held(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6))
        patch_held=internal_energy%fields(p)%values(1,:,:,:)
      endif
      if(linear_halos) then
      call fill_reactive_eb_tree_halo_3d(species,tree,p,ghost,temp,good,boundary_faces=boundary_faces, &
        inflow_state=inflow_state,inflow_temperature=inflow_temperature,inflow_mach=inflow_mach, &
        inflow_pressure=inflow_pressure,ghost_volume_fraction=vf, &
        halo_width=6,target_volume_fraction=extended_geometry(p)%volume_fraction,restrict_snapshot=.false., &
        ancestor_start=ancestor_start,ancestor_end=ancestor_end,ancestor_fraction=ancestor_fraction,frame=frame, &
        internal_energy=internal_energy,ghost_internal_energy=ghost_held, &
        ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end,coarse_geometry=extended_geometry, &
          separate_momentum=separate_momentum, &
        consistent_velocity=consistent_velocity,paired_pressure=paired_pressure,pressure_level=target_depth, &
      connected_directions=connected_directions,preserve_temperature_seed=use_total_energy)
      else
      call fill_reactive_eb_tree_halo_3d(species,tree,p,ghost,temp,good,boundary_faces=boundary_faces, &
        inflow_state=inflow_state,inflow_temperature=inflow_temperature,inflow_mach=inflow_mach, &
        inflow_pressure=inflow_pressure,ghost_volume_fraction=vf, &
        halo_width=6,target_volume_fraction=extended_geometry(p)%volume_fraction,restrict_snapshot=.false., &
        ancestor_start=ancestor_start,ancestor_end=ancestor_end,ancestor_fraction=ancestor_fraction,frame=frame, &
        internal_energy=internal_energy,ghost_internal_energy=ghost_held, &
        ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end)
      endif
      if(.not.good) return
      allocate(old%fields(p)%values(nv,n(1)+12,n(2)+12,n(3)+12))
      old%fields(p)%values=ghost
      ordinary=ghost
      if(.not.present(internal_energy)) then
      call ordinary_transport_energy(frame,tree%nodes(p)%geometry,[-5,-5,-5],vf>0,ordinary,good)
      endif
      if(.not.good) return
      level=tree%topology%nodes(p)%level+1
      faces=[origin(1,p)>0,origin(1,p)+n(1)<domain(1,level),origin(2,p)>0,origin(2,p)+n(2)<domain(2,level), &
        origin(3,p)>0,origin(3,p)+n(3)<domain(3,level)]
      if(present(inflow_mach)) faces=faces.or.(inflow_mach>=0)
      labels='outflow'
      if(present(boundary_faces)) then
        ! Patch edges within the level use explicit same-resolution donors.
        ! Boundary labels apply only where the patch touches a domain face.
        where(.not.faces) labels=boundary_faces
        faces=.true.
        do level=1,3
          if((labels(2*level-1)=='periodic').neqv.(labels(2*level)=='periodic')) labels(2*level-1:2*level)='outflow'
        enddo
      endif
      if(use_total_energy) then
      call compute_sdc_eb_transport_halo_3d(species,transport,ordinary,temp,extended_geometry(p),dt,flags, &
        rates%fields(p)%values,pending_flux(p),patch_halo,patch_theta,good,faces,labels,wall_temperature, &
        constant_coefficients=constant_coefficients)
      else
      call compute_sdc_eb_transport_halo_3d(species,transport,ordinary,temp,extended_geometry(p),dt,flags, &
        rates%fields(p)%values,pending_flux(p),patch_halo,patch_theta,good,faces,labels,wall_temperature, &
        constant_coefficients=constant_coefficients,internal_energy_density=ghost_held)
      endif
      if(.not.good) return
      theta=min(theta,patch_theta)
      if(wall_active) then
        allocate(body_halo(nv,n(1)+6,n(2)+6,n(3)+6),valid_body(nv,n(1),n(2),n(3)))
        call compute_eb_wall_transport_source_3d(species,transport,ordinary,temp,extended_geometry(p),[4,4,4],n+9, &
          present(embedded_wall_temperature),noslip,wall_t,body_halo,ignored_exchange,ignored_frequency,good, &
          constant_coefficients=constant_coefficients,internal_energy_density=ghost_held)
        if(.not.good) return
        call compute_eb_wall_transport_source_3d(species,transport,ordinary,temp,extended_geometry(p),[7,7,7],n+6, &
          present(embedded_wall_temperature),noslip,wall_t,valid_body,ignored_exchange,patch_frequency,good, &
          constant_coefficients=constant_coefficients,internal_energy_density=ghost_held)
        if(.not.good) return
        patch_halo(:,4:n(1)+9,4:n(2)+9,4:n(3)+9)=patch_halo(:,4:n(1)+9,4:n(2)+9,4:n(3)+9)+body_halo
        rates%fields(p)%values=rates%fields(p)%values+valid_body
        body%fields(p)%values=valid_body;frequency=max(frequency,patch_frequency)
        deallocate(body_halo,valid_body)
      endif
      if(.not.all(ieee_is_finite(patch_halo)).or..not.all(ieee_is_finite(rates%fields(p)%values))) return
      call move_alloc(patch_halo,halos%fields(p)%values)
      if(allocated(ghost_held)) deallocate(ghost_held)
      deallocate(ghost,temp,vf,ordinary)
    enddo
    call refresh_peers(good)
    if(.not.good) return
    if(wall_active) then
      ledger=tree
      do p=1,size(tree%nodes)
        ledger%nodes(p)%state=body%fields(p)%values
      enddo
      call integrate_reactive_eb_tree_3d(ledger,exchange,good)
      if(.not.good) return
      exchange=dt*exchange
      if(.not.all(ieee_is_finite(exchange))) return
    endif
    if(present(registers)) then
      first_depth=max(0,depth);last_depth=first_depth
      if(depth<0) last_depth=tree%topology%level_count()-1
      ! Non-MOL Diffusion.cpp ML StateRedist uses factor1. Stage flux
      ! acceptance (.5 for old/final) remains the outer SDC driver's job.
      do stage_depth=first_depth,last_depth
        call compute_state_redist_level_stage_3d(tree%topology,extended_geometry,old,halos,nv, &
          stage_depth,dt,rates,pending_registers,good,boundary_faces,target_volume_fraction,max_order, &
          register_scale=1._dp)
        if(.not.good) return
      enddo
    endif
    if(present(registers)) call move_alloc(pending_registers,registers)
    call move_alloc(rates%fields,rate%fields);call move_alloc(pending_flux,flux)
    call move_alloc(halos%fields,source_halos%fields);call move_alloc(old%fields,old_halos%fields)
    minimum_theta=theta
    if(wall_active) then
      embedded_exchange=exchange;embedded_frequency=frequency
    endif
    ok=.true.
  contains
    subroutine refresh_peers(valid)
      logical, intent(out) :: valid
      integer(int64) :: cell(3)
      integer :: p,q,i,j,k,a,level,index(3),m(3)
      logical :: inside
      character(len=32) :: faces(6)
      faces='outflow'
      if(present(boundary_faces)) faces=boundary_faces
      do p=1,size(tree%nodes)
        if(.not.allocated(halos%fields(p)%values)) cycle
        m=tree%topology%nodes(p)%shape;level=tree%topology%nodes(p)%level
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
          do q=1,size(tree%nodes)
            if(tree%topology%nodes(q)%level/=level) cycle
            if(any(cell<origin(:,q)).or.any(cell>=origin(:,q)+int(tree%topology%nodes(q)%shape,int64))) cycle
            index=int(cell-origin(:,q))+1
            halos%fields(p)%values(:,i,j,k)=rates%fields(q)%values(:,index(1),index(2),index(3))
            exit
          enddo
        enddo
        enddo
        enddo
      enddo
      valid=.true.
    end subroutine refresh_peers
  end subroutine compute_sdc_eb_tree_transport_halos_3d

  subroutine compute_sdc_eb_tree_transport_3d(species,transport,dt,viscosity,conduction, &
      diffusion,barodiffusion,tree,transport_rate,integrated_flux,minimum_theta,ok, &
      inflow_mach,inflow_temperature,inflow_pressure,boundary_faces,inflow_state,wall_temperature, &
      target_depth,ancestor_start,ancestor_end,ancestor_fraction, &
      embedded_wall_temperature,extended_geometry,embedded_wall_energy,embedded_frequency, &
      embedded_noslip,embedded_wall_impulse,frame, &
      constant_coefficients,internal_energy,ancestor_internal_energy_start,ancestor_internal_energy_end)
    real(dp), optional, intent(in) :: constant_coefficients(:)
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy,ancestor_internal_energy_start, &
      ancestor_internal_energy_end
    real(dp), allocatable :: ghost_held(:,:,:),patch_held(:,:,:)
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(in) :: dt
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(inout) :: transport_rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: integrated_flux(:)
    real(dp), intent(inout) :: minimum_theta
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: wall_temperature(6)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: target_depth
    type(eb_reactive_tree_3d), optional, intent(in) :: ancestor_start,ancestor_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    real(dp), optional, intent(in) :: embedded_wall_temperature
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry(:)
    real(dp), optional, intent(inout) :: embedded_wall_energy,embedded_frequency,embedded_wall_impulse(3)
    logical, optional, intent(in) :: embedded_noslip
    logical :: noslip,wall_active
    real(dp) :: wall_values(reactive_nvar(size(species))),heat_frequency,heat_energy,wall_impulse(3),wall_temp
    type(amr_tree_state_3d) :: heat
    integer :: selected_depth
    type(amr_tree_state_3d) :: rates
    type(eb_integrated_flux_3d), allocatable :: flux(:)
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),ordinary(:,:,:,:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp) :: theta,stage_theta
    logical :: faces(6)
    character(len=16) :: patch_faces(6)
    integer :: nv,node,n(3),depth
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if(present(internal_energy)) then
      if(.not.internal_energy%is_valid(tree%topology)) return
      if(size(internal_energy%fields(1)%values,1)/=1) return
    endif
    if(present(constant_coefficients)) then
      if(.not.valid_constant_transport_coefficients(constant_coefficients)) return
    endif
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    noslip=.false.;wall_temp=0
    if(present(embedded_noslip)) noslip=embedded_noslip
    wall_active=noslip.or.present(embedded_wall_temperature)
    if(noslip.neqv.present(embedded_wall_impulse)) return
    if(noslip.and..not.viscosity) return
    if(wall_active.neqv.present(extended_geometry)) return
    if(wall_active.neqv.present(embedded_wall_energy)) return
    if(wall_active.neqv.present(embedded_frequency)) return
    if(present(embedded_wall_temperature)) then
      wall_temp=embedded_wall_temperature
      if(.not.conduction.or..not.ieee_is_finite(wall_temp).or.wall_temp<=0) return
    end if
    if(present(wall_temperature)) then
      if(.not.present(boundary_faces)) return
      if(.not.all(ieee_is_finite(wall_temperature)).or.any(wall_temperature<=0)) return
    end if
    if (.not.ieee_is_finite(dt).or.dt<=0._dp) return
    if (.not.valid_tree(tree,nv)) return
    selected_depth=-1
    if(present(target_depth)) then
      selected_depth=target_depth
      if(selected_depth<0.or.selected_depth>=tree%topology%level_count()) return
    end if
    if(present(ancestor_start).or.present(ancestor_end).or.present(ancestor_fraction)) then
      if(selected_depth<0) return
      if(.not.(present(ancestor_start).and.present(ancestor_end).and.present(ancestor_fraction))) return
    end if

    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if (.not.valid) return
    allocate(flux(size(tree%nodes)),rates%fields(size(tree%nodes)));theta=1._dp
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      allocate(rates%fields(node)%values(nv,n(1),n(2),n(3)),source=0._dp)
      if(selected_depth>=0.and.tree%topology%nodes(node)%level/=selected_depth) cycle
      allocate(ghost(nv,0:n(1)+1,0:n(2)+1,0:n(3)+1),gt(0:n(1)+1,0:n(2)+1,0:n(3)+1))
      if(present(internal_energy)) then
        allocate(ghost_held(0:n(1)+1,0:n(2)+1,0:n(3)+1))
        patch_held=internal_energy%fields(node)%values(1,:,:,:)
      endif
      call fill_reactive_eb_tree_halo_3d(species,tree,node,ghost,gt,valid, &
        inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,restrict_snapshot=.false., &
        ancestor_start=ancestor_start,ancestor_end=ancestor_end,ancestor_fraction=ancestor_fraction,frame=frame, &
        internal_energy=internal_energy,ghost_internal_energy=ghost_held, &
        ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end)
      if (.not.valid) return
      depth=tree%topology%nodes(node)%level+1
      ! Inactive physical domain faces remain closed to diffusion. Parent patch edges
      ! inside the domain are donor interfaces, even when no peer exists.
      faces=[origin(1,node)>0,origin(1,node)+n(1)<domain(1,depth), &
        origin(2,node)>0,origin(2,node)+n(2)<domain(2,depth), &
        origin(3,node)>0,origin(3,node)+n(3)<domain(3,depth)]
      if (present(inflow_mach)) faces=faces.or.(inflow_mach>=0)
      patch_faces='outflow'
      if(present(boundary_faces)) then
        where(.not.faces.and.(boundary_faces=='slip_wall'.or.boundary_faces=='no_slip_wall'.or. &
              boundary_faces=='isothermal_wall')) &
          patch_faces=boundary_faces
        faces=.true.
      end if
      ordinary=tree%nodes(node)%state
      if(.not.present(internal_energy)) then
      call ordinary_transport_energy(frame,tree%nodes(node)%geometry,[1,1,1], &
        tree%nodes(node)%geometry%volume_fraction>0,ordinary,valid)
      endif
      if(.not.valid) return
      if(.not.present(internal_energy)) then
      call ordinary_transport_energy(frame,tree%nodes(node)%geometry,[0,0,0],gt>0,ghost,valid)
      endif
      if(.not.valid) return
      call compute_sdc_eb_transport_3d(species,transport,ordinary,tree%nodes(node)%temperature, &
        tree%nodes(node)%geometry,dt,[viscosity,conduction,diffusion,barodiffusion],rates%fields(node)%values, &
        flux(node),stage_theta,valid,ghost,gt,faces,physical_faces=patch_faces,wall_temperature=wall_temperature, &
        constant_coefficients=constant_coefficients,internal_energy_density=patch_held, &
        ghost_internal_energy_density=ghost_held)
      if (.not.valid) return
      theta=min(theta,stage_theta)
      if(allocated(ghost_held)) deallocate(ghost_held)
      deallocate(ghost,gt)
    end do
    if(.not.ieee_is_finite(theta)) return
    if(wall_active) then
      call compute_eb_tree_wall_transport_3d(species,transport,tree,extended_geometry, &
        present(embedded_wall_temperature),noslip,wall_temp,heat,wall_values,heat_frequency,valid, &
        boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure,target_depth, &
        ancestor_start,ancestor_end,ancestor_fraction,frame, &
        constant_coefficients=constant_coefficients,internal_energy=internal_energy, &
        ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end)
      if(.not.valid) return
      heat_energy=dt*wall_values(5);wall_impulse=dt*wall_values(2:4)
      if(.not.all(ieee_is_finite([heat_energy,wall_impulse]))) return
      do node=1,size(tree%nodes)
        rates%fields(node)%values=rates%fields(node)%values+heat%fields(node)%values
        if(.not.all(ieee_is_finite(rates%fields(node)%values))) return
      end do
    end if
    transport_rate=rates;call move_alloc(flux,integrated_flux)
    minimum_theta=theta;ok=.true.
    if(wall_active) then
      embedded_wall_energy=heat_energy;embedded_frequency=heat_frequency
      if(noslip) embedded_wall_impulse=wall_impulse
    end if
  end subroutine compute_sdc_eb_tree_transport_3d

  subroutine advance_reactive_eb_tree_transport_3d(species,transport,dt,viscosity,conduction, &
      diffusion,barodiffusion,tree,boundary_change,minimum_theta,ok, &
      inflow_mach,inflow_temperature,inflow_pressure,boundary_faces,inflow_state,wall_temperature, &
      constant_coefficients)
    real(dp), optional, intent(in) :: constant_coefficients(:)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(in) :: dt
    type(eb_reactive_tree_3d), intent(inout) :: tree
    real(dp), intent(inout) :: boundary_change(:),minimum_theta
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: wall_temperature(6)
    logical, intent(out) :: ok
    type(eb_reactive_tree_3d) :: snapshot,candidate
    type(eb_integrated_flux_3d), allocatable :: flux(:)
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),boundary(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp) :: theta,stage_theta
    logical :: faces(6)
    character(len=16) :: patch_faces(6)
    integer :: nv,node,n(3),depth,c
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if(present(constant_coefficients)) then
      if(.not.valid_constant_transport_coefficients(constant_coefficients)) return
    endif
    if(present(wall_temperature)) then
      if(.not.present(boundary_faces)) return
      if(.not.all(ieee_is_finite(wall_temperature)).or.any(wall_temperature<=0)) return
    end if
    if (.not.ieee_is_finite(dt).or.dt<=0._dp) return
    if (size(boundary_change)/=nv) return
    if (.not.valid_tree(tree,nv)) return
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if (.not.valid) return
    snapshot=tree
    call average_down_reactive_eb_tree_3d(species,snapshot,valid)
    if (.not.valid) return
    candidate=snapshot
    allocate(flux(size(tree%nodes)),boundary(nv));theta=1._dp
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      allocate(ghost(nv,0:n(1)+1,0:n(2)+1,0:n(3)+1),gt(0:n(1)+1,0:n(2)+1,0:n(3)+1))
      call fill_reactive_eb_tree_halo_3d(species,snapshot,node,ghost,gt,valid, &
        inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state)
      if (.not.valid) return
      depth=tree%topology%nodes(node)%level+1
      ! Inactive physical domain faces remain closed to diffusion. Parent patch edges
      ! inside the domain are donor interfaces, even when no peer exists.
      faces=[origin(1,node)>0,origin(1,node)+n(1)<domain(1,depth), &
        origin(2,node)>0,origin(2,node)+n(2)<domain(2,depth), &
        origin(3,node)>0,origin(3,node)+n(3)<domain(3,depth)]
      if (present(inflow_mach)) faces=faces.or.(inflow_mach>=0)
      patch_faces='outflow'
      if(present(boundary_faces)) then
        where(.not.faces.and.(boundary_faces=='slip_wall'.or.boundary_faces=='no_slip_wall'.or. &
              boundary_faces=='isothermal_wall')) &
          patch_faces=boundary_faces
        faces=.true.
      end if
      call advance_reactive_eb_transport_patch_euler_3d(species,transport,snapshot%nodes(node)%geometry, &
        dt,viscosity,conduction,diffusion,barodiffusion,candidate%nodes(node)%state, &
        candidate%nodes(node)%temperature,flux(node),stage_theta,valid,ghost,gt,faces, &
            physical_faces=patch_faces,wall_temperature=wall_temperature, &
        constant_coefficients=constant_coefficients)
      if (.not.valid) return
      theta=min(theta,stage_theta)
      deallocate(ghost,gt)
    end do
    call synchronize_reactive_eb_tree_flux_3d(species,candidate,flux,valid,boundary_faces=boundary_faces)
    if (.not.valid) return
    n=tree%topology%nodes(1)%shape
    associate(g=>tree%nodes(1)%geometry,f=>flux(1))
    do c=1,nv
      boundary(c)=g%dy*g%dz*sum(g%x_face_fraction(0,:,:)*f%x(c,0,:,:)- &
        g%x_face_fraction(n(1),:,:)*f%x(c,n(1),:,:))+ &
        g%dx*g%dz*sum(g%y_face_fraction(:,0,:)*f%y(c,:,0,:)- &
        g%y_face_fraction(:,n(2),:)*f%y(c,:,n(2),:))+ &
        g%dx*g%dy*sum(g%z_face_fraction(:,:,0)*f%z(c,:,:,0)- &
        g%z_face_fraction(:,:,n(3))*f%z(c,:,:,n(3)))
    end do
    end associate
    if (.not.all(ieee_is_finite(boundary))) return
    tree=candidate;boundary_change=boundary;minimum_theta=theta;ok=.true.
  end subroutine advance_reactive_eb_tree_transport_3d

  ! Serial common-interval hydro with PCM donors and outflow/Mach inlet edges.
  ! Supplying extended geometry selects PLM CTU; omission retains PCM Euler.
  ! Every node reads the same averaged snapshot; all outputs commit together.
  ! Provisional source-aware hydro: no restriction, redistribution or reflux.
  ! Rates are per fluid volume; fluxes are time-integrated face densities.
  subroutine compute_mol_eb_tree_hydro_3d(species,tree,extended_geometry,dt,order,hydro_rate,flux,ok, &
      boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
      target_depth,ancestor_start,ancestor_end,ancestor_fraction,lax_friedrichs,frame,spatial_rate_halos,old_state_halos, &
      internal_energy_rate,internal_energy_flux,internal_energy_spatial_rate_halos, &
      internal_energy_state,internal_energy_old_state_halos, &
      ancestor_internal_energy_start,ancestor_internal_energy_end, &
      separate_momentum,consistent_velocity,paired_pressure,connected_directions,diagnostics)
    type(eb_mol_diagnostics_3d), allocatable, optional, intent(inout) :: diagnostics(:)
    type(eb_mol_diagnostics_3d), allocatable :: pending_diagnostics(:),patch_diagnostics
    logical :: preserve_seed
    logical, optional, intent(in) :: separate_momentum,consistent_velocity,paired_pressure,connected_directions
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy_rate,internal_energy_spatial_rate_halos
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy_old_state_halos
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy_state,ancestor_internal_energy_start, &
      ancestor_internal_energy_end
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: internal_energy_flux(:)
    type(amr_tree_state_3d), optional, intent(inout) :: spatial_rate_halos,old_state_halos
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(eb_geometry_3d), intent(in) :: extended_geometry(:)
    real(dp), intent(in) :: dt
    integer, intent(in) :: order
    type(amr_tree_state_3d), intent(inout) :: hydro_rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: lax_friedrichs
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    integer, optional, intent(in) :: target_depth
    type(eb_reactive_tree_3d), optional, intent(in) :: ancestor_start,ancestor_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    type(amr_tree_state_3d) :: rates,alpha,pending_halos,pending_old,auxiliary_rates,pending_aux_halos,pending_aux_old
    type(eb_integrated_flux_3d), allocatable :: candidate_flux(:),raw_flux(:),auxiliary_flux(:),auxiliary_raw(:)
    type(eb_geometry_3d) :: center
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),vf(:,:,:),patch_halo(:,:,:,:),patch_aux_halo(:,:,:),ghost_aux(:,:,:)
    real(dp) :: boundary(reactive_nvar(size(species))),wall(reactive_nvar(size(species)))
    integer :: nv,node,n(3),depth,halo_width
    logical :: valid,use_lax,use_total_energy
    ! The ordinary nonrotating path follows pinned pc_ctoprim for both face
    ! and wall pressure. Explicit experimental halo policies retain their
    ! held-energy pressure definition, as does rotating cancellation recovery.
    use_total_energy=.not.present(frame)
    if(present(separate_momentum)) then
      if(separate_momentum) use_total_energy=.false.
    endif
    preserve_seed=.not.present(frame)
    if(present(separate_momentum)) preserve_seed=preserve_seed.and..not.separate_momentum
    ok=.false.;nv=reactive_nvar(size(species));depth=-1
    if(present(consistent_velocity)) then
      if(consistent_velocity) then
        if(.not.present(separate_momentum)) return
        if(.not.separate_momentum) return
      endif
    endif
    if(present(connected_directions)) then
      if(connected_directions) then
        if(.not.present(separate_momentum)) return
        if(.not.separate_momentum) return
      endif
    endif
    if(present(separate_momentum)) then
      if(separate_momentum.and..not.present(internal_energy_old_state_halos)) return
    endif
    if(present(internal_energy_rate).neqv.present(internal_energy_flux)) return
    if(present(internal_energy_spatial_rate_halos).and..not.present(internal_energy_rate)) return
    if(present(internal_energy_state).neqv.present(internal_energy_old_state_halos)) return
    if(present(internal_energy_state).and..not.present(internal_energy_rate)) return
    if(present(ancestor_internal_energy_start).neqv.present(ancestor_internal_energy_end)) return
    if(present(ancestor_internal_energy_start).and..not.present(internal_energy_state)) return
    if(present(internal_energy_state)) then
      if(present(ancestor_fraction).neqv.present(ancestor_internal_energy_start)) return
    endif
    use_lax=.false.
    if(present(lax_friedrichs)) use_lax=lax_friedrichs
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    if(.not.valid_tree(tree,nv).or..not.ieee_is_finite(dt)) return
    if(dt<=0.or.(order/=1.and.order/=2)) return
    if(size(extended_geometry)/=size(tree%nodes)) return
    if(present(target_depth)) then
      depth=target_depth
      if(depth<0.or.depth>=tree%topology%level_count()) return
    end if
    if(present(ancestor_start).or.present(ancestor_end).or.present(ancestor_fraction)) then
      if(depth<0) return
    end if
    halo_width=3
    if(present(spatial_rate_halos).or.present(old_state_halos).or. &
      present(internal_energy_spatial_rate_halos).or.present(internal_energy_old_state_halos)) halo_width=6
    allocate(rates%fields(size(tree%nodes)),candidate_flux(size(tree%nodes)))
    if(present(diagnostics)) allocate(pending_diagnostics(size(tree%nodes)),patch_diagnostics)
    if(present(internal_energy_rate)) then
      allocate(auxiliary_rates%fields(size(tree%nodes)),auxiliary_flux(size(tree%nodes)))
      if(use_lax) allocate(auxiliary_raw(size(tree%nodes)))
    endif
    if(halo_width==6) allocate(pending_halos%fields(size(tree%nodes)),pending_old%fields(size(tree%nodes)))
    if(present(internal_energy_spatial_rate_halos)) allocate(pending_aux_halos%fields(size(tree%nodes)))
    if(present(internal_energy_old_state_halos)) allocate(pending_aux_old%fields(size(tree%nodes)))
    if(use_lax) allocate(alpha%fields(size(tree%nodes)),raw_flux(size(tree%nodes)))
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      if(any([extended_geometry(node)%nx,extended_geometry(node)%ny,extended_geometry(node)%nz]/=n+2*halo_width)) return
      call extract_eb_geometry_3d(extended_geometry(node),[1,1,1]+halo_width,n+halo_width,center,valid)
      if(.not.valid) return
      if(.not.matching_eb_metrics_3d(center,tree%nodes(node)%geometry)) return
      allocate(rates%fields(node)%values(nv,n(1),n(2),n(3)),source=0._dp)
      if(present(internal_energy_rate)) &
        allocate(auxiliary_rates%fields(node)%values(1,n(1),n(2),n(3)),source=0._dp)
      if(use_lax) then
        allocate(alpha%fields(node)%values(nv-5,n(1),n(2),n(3)),source=1._dp)
      end if
      if(depth>=0.and.tree%topology%nodes(node)%level/=depth) cycle
      allocate(ghost(nv,1-halo_width:n(1)+halo_width,1-halo_width:n(2)+halo_width,1-halo_width:n(3)+halo_width), &
        gt(1-halo_width:n(1)+halo_width,1-halo_width:n(2)+halo_width,1-halo_width:n(3)+halo_width), &
        vf(1-halo_width:n(1)+halo_width,1-halo_width:n(2)+halo_width,1-halo_width:n(3)+halo_width))
      if(present(internal_energy_old_state_halos)) then
        allocate(ghost_aux(1-halo_width:n(1)+halo_width,1-halo_width:n(2)+halo_width,1-halo_width:n(3)+halo_width))
        call fill_reactive_eb_tree_halo_3d(species,tree,node,ghost,gt,valid, &
          boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
          inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,ghost_volume_fraction=vf, &
          halo_width=halo_width,target_volume_fraction=extended_geometry(node)%volume_fraction,restrict_snapshot=.false., &
          ancestor_start=ancestor_start,ancestor_end=ancestor_end,ancestor_fraction=ancestor_fraction,frame=frame, &
          internal_energy=internal_energy_state,ghost_internal_energy=ghost_aux, &
          ancestor_internal_energy_start=ancestor_internal_energy_start, &
          ancestor_internal_energy_end=ancestor_internal_energy_end,coarse_geometry=extended_geometry, &
          separate_momentum=separate_momentum, &
        consistent_velocity=consistent_velocity,paired_pressure=paired_pressure,pressure_level=target_depth, &
      connected_directions=connected_directions,preserve_temperature_seed=preserve_seed)
      else
        call fill_reactive_eb_tree_halo_3d(species,tree,node,ghost,gt,valid, &
          boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
          inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,ghost_volume_fraction=vf, &
          halo_width=halo_width,target_volume_fraction=extended_geometry(node)%volume_fraction,restrict_snapshot=.false., &
          ancestor_start=ancestor_start,ancestor_end=ancestor_end,ancestor_fraction=ancestor_fraction,frame=frame)
      endif
      if(.not.valid) return
      if(present(internal_energy_old_state_halos)) then
        allocate(pending_aux_old%fields(node)%values(1,n(1)+12,n(2)+12,n(3)+12))
        pending_aux_old%fields(node)%values(1,:,:,:)=ghost_aux
      endif
      if(halo_width==6) then
        if(use_lax) then
          if(present(internal_energy_rate)) then
            call compute_mol_eb_hydro_3d(species,extended_geometry(node),ghost,gt,dt,order, &
              rates%fields(node)%values,candidate_flux(node),boundary,wall,valid,lax_friedrichs,frame, &
              raw_flux(node),alpha%fields(node)%values,.true.,patch_halo, &
              internal_energy_rate=auxiliary_rates%fields(node)%values(1,:,:,:), &
              internal_energy_integrated_flux=auxiliary_flux(node),internal_energy_raw_face_flux=auxiliary_raw(node), &
              internal_energy_spatial_rate_halo=patch_aux_halo,internal_energy_density=ghost_aux, &
              diagnostics=patch_diagnostics,primitive_from_total_energy=use_total_energy)
          else
            call compute_mol_eb_hydro_3d(species,extended_geometry(node),ghost,gt,dt,order, &
              rates%fields(node)%values,candidate_flux(node),boundary,wall,valid,lax_friedrichs,frame, &
              raw_flux(node),alpha%fields(node)%values,.true.,patch_halo, &
              diagnostics=patch_diagnostics,primitive_from_total_energy=use_total_energy)
          endif
        else
          if(present(internal_energy_rate)) then
            call compute_mol_eb_hydro_3d(species,extended_geometry(node),ghost,gt,dt,order, &
              rates%fields(node)%values,candidate_flux(node),boundary,wall,valid,lax_friedrichs,frame, &
              spatial_rate_halo=patch_halo,internal_energy_rate=auxiliary_rates%fields(node)%values(1,:,:,:), &
              internal_energy_integrated_flux=auxiliary_flux(node), &
              internal_energy_spatial_rate_halo=patch_aux_halo,internal_energy_density=ghost_aux, &
              diagnostics=patch_diagnostics,primitive_from_total_energy=use_total_energy)
          else
            call compute_mol_eb_hydro_3d(species,extended_geometry(node),ghost,gt,dt,order, &
              rates%fields(node)%values,candidate_flux(node),boundary,wall,valid,lax_friedrichs,frame, &
      spatial_rate_halo=patch_halo, &
              diagnostics=patch_diagnostics,primitive_from_total_energy=use_total_energy)
          endif
        endif
        if(.not.valid) return
        call move_alloc(patch_halo,pending_halos%fields(node)%values)
        if(present(internal_energy_spatial_rate_halos)) then
          allocate(pending_aux_halos%fields(node)%values(1,size(patch_aux_halo,1), &
            size(patch_aux_halo,2),size(patch_aux_halo,3)))
          pending_aux_halos%fields(node)%values(1,:,:,:)=patch_aux_halo
        endif
        if(allocated(patch_aux_halo)) deallocate(patch_aux_halo)
        allocate(pending_old%fields(node)%values(nv,n(1)+12,n(2)+12,n(3)+12))
        pending_old%fields(node)%values=ghost
      else
      if(use_lax) then
        if(present(internal_energy_rate)) then
          call compute_mol_eb_hydro_3d(species,extended_geometry(node),ghost,gt,dt,order, &
            rates%fields(node)%values,candidate_flux(node),boundary,wall,valid,lax_friedrichs,frame, &
            raw_flux(node),alpha%fields(node)%values,.true., &
            internal_energy_rate=auxiliary_rates%fields(node)%values(1,:,:,:), &
            internal_energy_integrated_flux=auxiliary_flux(node),internal_energy_raw_face_flux=auxiliary_raw(node), &
            internal_energy_density=ghost_aux, &
              diagnostics=patch_diagnostics,primitive_from_total_energy=use_total_energy)
        else
          call compute_mol_eb_hydro_3d(species,extended_geometry(node),ghost,gt,dt,order, &
            rates%fields(node)%values,candidate_flux(node),boundary,wall,valid,lax_friedrichs,frame, &
            raw_flux(node),alpha%fields(node)%values,.true., &
              diagnostics=patch_diagnostics,primitive_from_total_energy=use_total_energy)
        endif
      else
        if(present(internal_energy_rate)) then
          call compute_mol_eb_hydro_3d(species,extended_geometry(node),ghost,gt,dt,order, &
            rates%fields(node)%values,candidate_flux(node),boundary,wall,valid,lax_friedrichs,frame, &
            internal_energy_rate=auxiliary_rates%fields(node)%values(1,:,:,:), &
            internal_energy_integrated_flux=auxiliary_flux(node),internal_energy_density=ghost_aux, &
              diagnostics=patch_diagnostics,primitive_from_total_energy=use_total_energy)
        else
          call compute_mol_eb_hydro_3d(species,extended_geometry(node),ghost,gt,dt,order, &
            rates%fields(node)%values,candidate_flux(node),boundary,wall,valid,lax_friedrichs,frame, &
              diagnostics=patch_diagnostics,primitive_from_total_energy=use_total_energy)
        endif
      end if
      endif
      if(.not.valid) return
      if(present(diagnostics)) pending_diagnostics(node)=patch_diagnostics
      if(allocated(ghost_aux)) deallocate(ghost_aux)
      deallocate(ghost,gt,vf)
    end do
    if(use_lax) then
      call limit_boundary_faces(valid)
      if(.not.valid) return
    end if
    if(halo_width==6) then
      ! Cross-patch centroid limiting may have adjusted valid rates. Keep
      ! the exported valid slice identical to the published spatial rate.
      do node=1,size(tree%nodes)
        if(.not.allocated(pending_halos%fields(node)%values)) cycle
        n=tree%topology%nodes(node)%shape
        pending_halos%fields(node)%values(:,7:n(1)+6,7:n(2)+6,7:n(3)+6)=rates%fields(node)%values
        if(present(internal_energy_spatial_rate_halos)) &
          pending_aux_halos%fields(node)%values(1,7:n(1)+6,7:n(2)+6,7:n(3)+6)= &
          auxiliary_rates%fields(node)%values(1,:,:,:)
      enddo
      call refresh_peer_source_halos(valid)
      if(.not.valid) return
    endif
    call move_alloc(rates%fields,hydro_rate%fields);call move_alloc(candidate_flux,flux)
    if(present(internal_energy_rate)) then
      call move_alloc(auxiliary_rates%fields,internal_energy_rate%fields)
      call move_alloc(auxiliary_flux,internal_energy_flux)
    endif
    if(present(spatial_rate_halos)) call move_alloc(pending_halos%fields,spatial_rate_halos%fields)
    if(present(old_state_halos)) call move_alloc(pending_old%fields,old_state_halos%fields)
    if(present(internal_energy_spatial_rate_halos)) &
      call move_alloc(pending_aux_halos%fields,internal_energy_spatial_rate_halos%fields)
    if(present(internal_energy_old_state_halos)) &
      call move_alloc(pending_aux_old%fields,internal_energy_old_state_halos%fields)
    if(present(diagnostics)) call move_alloc(pending_diagnostics,diagnostics)
    ok=.true.
  contains
    ! After centroid-flux limiting, peer ghost rates must equal the peer's
    ! accepted valid spatial rate. Coarse-interpolated ghosts retain the RHS
    ! evaluated locally at this level; no coarse-rate interpolation is used.
    subroutine refresh_peer_source_halos(success)
      logical, intent(out) :: success
      integer(int64), allocatable :: origin(:,:),domain(:,:)
      integer, allocatable :: ratios(:)
      integer(int64) :: cell(3)
      integer :: p,q,level,i,j,k,index(3),a,m(3)
      character(len=32) :: faces(6)
      logical :: inside
      faces='outflow'
      if(present(boundary_faces)) faces=boundary_faces
      call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,success)
      if(.not.success) return
      do p=1,size(tree%nodes)
        if(.not.allocated(pending_halos%fields(p)%values)) cycle
        level=tree%topology%nodes(p)%level;m=tree%topology%nodes(p)%shape
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
          do q=1,size(tree%nodes)
            if(tree%topology%nodes(q)%level/=level) cycle
            if(any(cell<origin(:,q)).or.any(cell>=origin(:,q)+int(tree%topology%nodes(q)%shape,int64))) cycle
            index=int(cell-origin(:,q))+1
            pending_halos%fields(p)%values(:,i,j,k)=rates%fields(q)%values(:,index(1),index(2),index(3))
            if(present(internal_energy_spatial_rate_halos)) &
              pending_aux_halos%fields(p)%values(1,i,j,k)= &
              auxiliary_rates%fields(q)%values(1,index(1),index(2),index(3))
            exit
          enddo
        enddo
        enddo
        enddo
      enddo
    end subroutine refresh_peer_source_halos

    subroutine limit_boundary_faces(success)
      logical, intent(out) :: success
      integer(int64), allocatable :: origin(:,:),domain(:,:)
      integer, allocatable :: ratios(:)
      integer(int64) :: global_cell(3),coarse_cell(3)
      integer :: p,q,axis,side,i,j,k,m(3),level,face(3),inside(3),neighbor(3),peer,s
      real(dp) :: neighbor_alpha(nv-5),delta(nv),delta_aux,theta,area,volume,sign_face
      character(len=16) :: faces(6)
      success=.false.;faces='outflow'
      if(present(boundary_faces)) faces=boundary_faces
      call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,success)
      if(.not.success) return
      do p=1,size(tree%nodes)
        if(depth>=0.and.tree%topology%nodes(p)%level/=depth) cycle
        m=tree%topology%nodes(p)%shape;level=tree%topology%nodes(p)%level+1
        do axis=1,3
        do side=0,1
          do k=1,m(3)
          do j=1,m(2)
          do i=1,m(1)
            face=[i,j,k]
            if(side==0.and.face(axis)/=1) cycle
            if(side==1.and.face(axis)/=m(axis)) cycle
            if(side==0) face(axis)=0
            inside=face
            if(side==0) inside(axis)=1
            global_cell=origin(:,p)+int(inside-1,int64)
            global_cell(axis)=global_cell(axis)+int(2*side-1,int64)
            if(faces(2*axis-1)=='periodic'.and.faces(2*axis)=='periodic') &
              global_cell(axis)=modulo(global_cell(axis),domain(axis,level))
            peer=0
            do q=1,size(tree%nodes)
              if(tree%topology%nodes(q)%level/=level-1) cycle
              if(all(global_cell>=origin(:,q)).and. &
                all(global_cell<origin(:,q)+int(tree%topology%nodes(q)%shape,int64))) then
                peer=q;exit
              end if
            end do
            if(peer==0.and.level>1) then
              if(all(global_cell>=0_int64).and.all(global_cell<domain(:,level))) then
                coarse_cell=global_cell/int(ratios(level),int64)
                do q=1,size(tree%nodes)
                  if(tree%topology%nodes(q)%level/=level-2) cycle
                  if(all(coarse_cell>=origin(:,q)).and. &
                    all(coarse_cell<origin(:,q)+int(tree%topology%nodes(q)%shape,int64))) then
                    peer=q;exit
                  end if
                end do
              end if
            end if
            neighbor_alpha=1._dp
            if(peer>0) then
              if(tree%topology%nodes(peer)%level==level-1) then
                neighbor=int(global_cell-origin(:,peer))+1
              else
                neighbor=int(coarse_cell-origin(:,peer))+1
              end if
              neighbor_alpha=alpha%fields(peer)%values(:,neighbor(1),neighbor(2),neighbor(3))
            end if
            select case(axis)
            case(1)
              delta=candidate_flux(p)%x(:,face(1),face(2),face(3))-raw_flux(p)%x(:,face(1),face(2),face(3))
              if(present(internal_energy_rate)) delta_aux=auxiliary_flux(p)%x(1,face(1),face(2),face(3))- &
                auxiliary_raw(p)%x(1,face(1),face(2),face(3))
              area=tree%nodes(p)%geometry%dy*tree%nodes(p)%geometry%dz* &
                tree%nodes(p)%geometry%x_face_fraction(face(1),face(2),face(3))
            case(2)
              delta=candidate_flux(p)%y(:,face(1),face(2),face(3))-raw_flux(p)%y(:,face(1),face(2),face(3))
              if(present(internal_energy_rate)) delta_aux=auxiliary_flux(p)%y(1,face(1),face(2),face(3))- &
                auxiliary_raw(p)%y(1,face(1),face(2),face(3))
              area=tree%nodes(p)%geometry%dx*tree%nodes(p)%geometry%dz* &
                tree%nodes(p)%geometry%y_face_fraction(face(1),face(2),face(3))
            case(3)
              delta=candidate_flux(p)%z(:,face(1),face(2),face(3))-raw_flux(p)%z(:,face(1),face(2),face(3))
              if(present(internal_energy_rate)) delta_aux=auxiliary_flux(p)%z(1,face(1),face(2),face(3))- &
                auxiliary_raw(p)%z(1,face(1),face(2),face(3))
              area=tree%nodes(p)%geometry%dx*tree%nodes(p)%geometry%dy* &
                tree%nodes(p)%geometry%z_face_fraction(face(1),face(2),face(3))
            end select
            if(all(delta(6:)==0._dp)) cycle
            theta=1._dp
            do s=1,nv-5
              if(delta(s+5)>0) then
                if(side==1) theta=min(theta,alpha%fields(p)%values(s,inside(1),inside(2),inside(3)))
                if(side==0) theta=min(theta,neighbor_alpha(s))
              else if(delta(s+5)<0) then
                if(side==0) theta=min(theta,alpha%fields(p)%values(s,inside(1),inside(2),inside(3)))
                if(side==1) theta=min(theta,neighbor_alpha(s))
              end if
            end do
            if(theta>=1._dp) cycle
            volume=tree%nodes(p)%geometry%volume_fraction(inside(1),inside(2),inside(3))* &
              tree%nodes(p)%geometry%dx*tree%nodes(p)%geometry%dy*tree%nodes(p)%geometry%dz
            if(.not.ieee_is_finite(volume).or.volume<=0._dp) then
              success=.false.;return
            end if
            sign_face=merge(1._dp,-1._dp,side==0)
            rates%fields(p)%values(:,inside(1),inside(2),inside(3))= &
              rates%fields(p)%values(:,inside(1),inside(2),inside(3))+ &
              sign_face*area*(theta-1._dp)*delta/(dt*volume)
            if(present(internal_energy_rate)) auxiliary_rates%fields(p)%values(1,inside(1),inside(2),inside(3))= &
              auxiliary_rates%fields(p)%values(1,inside(1),inside(2),inside(3))+ &
              sign_face*area*(theta-1._dp)*delta_aux/(dt*volume)
            select case(axis)
            case(1)
              candidate_flux(p)%x(:,face(1),face(2),face(3))= &
                raw_flux(p)%x(:,face(1),face(2),face(3))+theta*delta
              if(present(internal_energy_rate)) auxiliary_flux(p)%x(1,face(1),face(2),face(3))= &
                auxiliary_raw(p)%x(1,face(1),face(2),face(3))+theta*delta_aux
            case(2)
              candidate_flux(p)%y(:,face(1),face(2),face(3))= &
                raw_flux(p)%y(:,face(1),face(2),face(3))+theta*delta
              if(present(internal_energy_rate)) auxiliary_flux(p)%y(1,face(1),face(2),face(3))= &
                auxiliary_raw(p)%y(1,face(1),face(2),face(3))+theta*delta_aux
            case(3)
              candidate_flux(p)%z(:,face(1),face(2),face(3))= &
                raw_flux(p)%z(:,face(1),face(2),face(3))+theta*delta
              if(present(internal_energy_rate)) auxiliary_flux(p)%z(1,face(1),face(2),face(3))= &
                auxiliary_raw(p)%z(1,face(1),face(2),face(3))+theta*delta_aux
            end select
          end do
          end do
          end do
        end do
        end do
        if(.not.all(ieee_is_finite(rates%fields(p)%values))) then
          success=.false.;return
        end if
      end do
    end subroutine limit_boundary_faces
  end subroutine compute_mol_eb_tree_hydro_3d

  ! Same-level conservative correction for complete forward source budgets.
  ! Input faces are dt-integrated and must agree on shared/periodic seams.
  ! Coarse/fine reflux retains its separate composite donor requirement.
  subroutine apply_eb_tree_flux_correction_3d(tree,source,baseline,dt,hydro_rate,flux,ok,periodic,target_depth, &
      internal_energy_baseline,internal_energy_rate,internal_energy_flux)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: source
    type(eb_integrated_flux_3d), intent(in) :: baseline(:)
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: hydro_rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: periodic(3)
    integer, optional, intent(in) :: target_depth
    type(eb_integrated_flux_3d), optional, intent(in) :: internal_energy_baseline(:)
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy_rate
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: internal_energy_flux(:)
    type(amr_tree_state_3d) :: held_candidate
    type(eb_integrated_flux_3d), allocatable :: held_limited(:)
    type(eb_integrated_flux_3d) :: held_low,held_original,held_delta
    real(dp) :: held_increment
    type(amr_tree_state_3d) :: candidate,allowance
    type(eb_integrated_flux_3d), allocatable :: original(:),limited(:)
    type(eb_integrated_flux_3d) :: low,delta
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: halo(:,:,:,:),endpoint(:)
    integer(int64) :: global_cell(3)
    real(dp) :: volume,areas(3)
    integer :: nv,nc,p,q,depth,n(3),level,axis,side,i,j,k,c,cell(3),neighbor(3)
    logical :: valid,wrap(3)
    ok=.false.;depth=-1;wrap=.false.
    if(present(periodic)) wrap=periodic
    if(present(target_depth)) depth=target_depth
    if(.not.ieee_is_finite(dt).or.dt<=0.or.depth< -1) return
    if(.not.allocated(tree%nodes).or..not.allocated(hydro_rate%fields).or..not.allocated(source%fields)) return
    if(.not.allocated(flux)) return
    if(size(hydro_rate%fields)/=size(tree%nodes).or.size(source%fields)/=size(tree%nodes).or. &
      size(flux)/=size(tree%nodes).or.size(baseline)/=size(tree%nodes).or.size(tree%nodes)<1) return
    if(.not.allocated(tree%nodes(1)%state)) return
    nv=size(tree%nodes(1)%state,1);nc=nv-4
    if(nv<6.or..not.valid_tree(tree,nv).or.depth>=tree%topology%level_count()) return
    do p=1,size(tree%nodes)
      n=tree%topology%nodes(p)%shape
      if(.not.allocated(hydro_rate%fields(p)%values).or..not.allocated(source%fields(p)%values)) return
      if(any(shape(hydro_rate%fields(p)%values)/=[nv,n]).or.any(shape(source%fields(p)%values)/=[nv,n])) return
    enddo
    if(present(internal_energy_baseline).neqv.present(internal_energy_rate)) return
    if(present(internal_energy_rate).neqv.present(internal_energy_flux)) return
    if(present(internal_energy_rate)) then
      if(.not.internal_energy_rate%is_valid(tree%topology)) return
      if(.not.allocated(internal_energy_flux)) return
      if(size(internal_energy_baseline)/=size(tree%nodes).or.size(internal_energy_flux)/=size(tree%nodes)) return
      do p=1,size(tree%nodes)
        if(size(internal_energy_rate%fields(p)%values,1)/=1) return
        if(.not.all(ieee_is_finite(internal_energy_rate%fields(p)%values))) return
        if(depth>=0.and.tree%topology%nodes(p)%level/=depth) cycle
        if(.not.valid_eb_integrated_flux_3d(internal_energy_baseline(p),tree%nodes(p)%geometry,1).or. &
          .not.valid_eb_integrated_flux_3d(internal_energy_flux(p),tree%nodes(p)%geometry,1)) return
      enddo
      held_candidate=internal_energy_rate
      allocate(held_limited(size(tree%nodes)))
    endif
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if(.not.valid) return
    candidate=hydro_rate
    allocate(allowance%fields(size(tree%nodes)),original(size(tree%nodes)),limited(size(tree%nodes)),endpoint(nv))
    do p=1,size(tree%nodes)
      n=tree%topology%nodes(p)%shape
      allocate(allowance%fields(p)%values(nc,n(1),n(2),n(3)),source=1._dp)
      if(depth>=0.and.tree%topology%nodes(p)%level/=depth) cycle
      original(p)=flux(p);low=baseline(p)
      if(.not.allocated(original(p)%x).or..not.allocated(original(p)%y).or..not.allocated(original(p)%z)) return
      if(.not.allocated(low%x).or..not.allocated(low%y).or..not.allocated(low%z)) return
      original(p)%x=original(p)%x/dt;original(p)%y=original(p)%y/dt;original(p)%z=original(p)%z/dt
      low%x=low%x/dt;low%y=low%y/dt;low%z=low%z/dt
      limited(p)=original(p)
      call limit_reactive_eb_flux_correction_3d(tree%nodes(p)%geometry,tree%nodes(p)%state,low,limited(p),dt,valid, &
        source%fields(p)%values,allowance%fields(p)%values)
      if(.not.valid) return
    enddo
    do p=1,size(tree%nodes)
      if(depth>=0.and.tree%topology%nodes(p)%level/=depth) cycle
      n=tree%topology%nodes(p)%shape;level=tree%topology%nodes(p)%level+1
      allocate(halo(nc,0:n(1)+1,0:n(2)+1,0:n(3)+1));halo=1
      halo(:,1:n(1),1:n(2),1:n(3))=allowance%fields(p)%values
      do axis=1,3
      do side=0,1
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        cell=[i,j,k]
        if(side==0.and.cell(axis)/=1) cycle
        if(side==1.and.cell(axis)/=n(axis)) cycle
        cell(axis)=cell(axis)+2*side-1
        global_cell=origin(:,p)+int(cell-1,int64)
        if(wrap(axis)) global_cell(axis)=modulo(global_cell(axis),domain(axis,level))
        do q=1,size(tree%nodes)
          if(tree%topology%nodes(q)%level/=tree%topology%nodes(p)%level) cycle
          if(all(global_cell>=origin(:,q)).and.all(global_cell<origin(:,q)+int(tree%topology%nodes(q)%shape,int64))) then
            neighbor=int(global_cell-origin(:,q))+1
            halo(:,cell(1),cell(2),cell(3))=allowance%fields(q)%values(:,neighbor(1),neighbor(2),neighbor(3))
            exit
          endif
        enddo
      enddo
      enddo
      enddo
      enddo
      enddo
      low=baseline(p);low%x=low%x/dt;low%y=low%y/dt;low%z=low%z/dt;limited(p)=original(p)
      if(present(internal_energy_rate)) then
        held_low=internal_energy_baseline(p)
        held_low%x=held_low%x/dt;held_low%y=held_low%y/dt;held_low%z=held_low%z/dt
        held_original=internal_energy_flux(p)
        held_original%x=held_original%x/dt;held_original%y=held_original%y/dt;held_original%z=held_original%z/dt
        held_limited(p)=held_original
        call limit_reactive_eb_flux_correction_3d(tree%nodes(p)%geometry,tree%nodes(p)%state,low,limited(p),dt,valid, &
          source=source%fields(p)%values,peer_alpha=halo, &
          internal_energy_baseline=held_low,internal_energy_flux=held_limited(p))
      else
      call limit_reactive_eb_flux_correction_3d(tree%nodes(p)%geometry,tree%nodes(p)%state,low,limited(p),dt,valid, &
        source=source%fields(p)%values,peer_alpha=halo)
      endif
      if(.not.valid) return
      deallocate(halo)
      where(limited(p)%x==original(p)%x)
        limited(p)%x=flux(p)%x
      elsewhere
        limited(p)%x=dt*limited(p)%x
      endwhere
      where(limited(p)%y==original(p)%y)
        limited(p)%y=flux(p)%y
      elsewhere
        limited(p)%y=dt*limited(p)%y
      endwhere
      where(limited(p)%z==original(p)%z)
        limited(p)%z=flux(p)%z
      elsewhere
        limited(p)%z=dt*limited(p)%z
      endwhere
      if(present(internal_energy_rate)) then
        where(held_limited(p)%x==held_original%x)
          held_limited(p)%x=internal_energy_flux(p)%x
        elsewhere
          held_limited(p)%x=dt*held_limited(p)%x
        endwhere
        where(held_limited(p)%y==held_original%y)
          held_limited(p)%y=internal_energy_flux(p)%y
        elsewhere
          held_limited(p)%y=dt*held_limited(p)%y
        endwhere
        where(held_limited(p)%z==held_original%z)
          held_limited(p)%z=internal_energy_flux(p)%z
        elsewhere
          held_limited(p)%z=dt*held_limited(p)%z
        endwhere
        if(.not.valid_eb_integrated_flux_3d(held_limited(p),tree%nodes(p)%geometry,1)) then
          return
        endif
        held_delta=held_limited(p)
        held_delta%x=held_delta%x-internal_energy_flux(p)%x
        held_delta%y=held_delta%y-internal_energy_flux(p)%y
        held_delta%z=held_delta%z-internal_energy_flux(p)%z
      endif
      delta=limited(p);delta%x=delta%x-flux(p)%x;delta%y=delta%y-flux(p)%y;delta%z=delta%z-flux(p)%z
      associate(g=>tree%nodes(p)%geometry)
      areas=[g%dy*g%dz,g%dx*g%dz,g%dx*g%dy]
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(g%volume_fraction(i,j,k)<=0) cycle
        volume=dt*g%dx*g%dy*g%dz*g%volume_fraction(i,j,k)
        endpoint=( &
          areas(1)*(g%x_face_fraction(i-1,j,k)*delta%x(:,i-1,j,k)-g%x_face_fraction(i,j,k)*delta%x(:,i,j,k))+ &
          areas(2)*(g%y_face_fraction(i,j-1,k)*delta%y(:,i,j-1,k)-g%y_face_fraction(i,j,k)*delta%y(:,i,j,k))+ &
          areas(3)*(g%z_face_fraction(i,j,k-1)*delta%z(:,i,j,k-1)-g%z_face_fraction(i,j,k)*delta%z(:,i,j,k)))/volume
        if(any(endpoint/=0)) candidate%fields(p)%values(:,i,j,k)=hydro_rate%fields(p)%values(:,i,j,k)+endpoint
        ! Slip walls carry no density or species flux. Recompute changed
        ! components directly: high divergence plus its opposite correction
        ! can leave a negative rounding residual when a species flux is zero.
        ! Retain the original rate for unchanged components and wall terms.
        do c=1,nv
          if(c>=2.and.c<=5) cycle
          if(endpoint(c)==0) cycle
          candidate%fields(p)%values(c,i,j,k)=( &
            areas(1)*(g%x_face_fraction(i-1,j,k)*limited(p)%x(c,i-1,j,k)- &
              g%x_face_fraction(i,j,k)*limited(p)%x(c,i,j,k))+ &
            areas(2)*(g%y_face_fraction(i,j-1,k)*limited(p)%y(c,i,j-1,k)- &
              g%y_face_fraction(i,j,k)*limited(p)%y(c,i,j,k))+ &
            areas(3)*(g%z_face_fraction(i,j,k-1)*limited(p)%z(c,i,j,k-1)- &
              g%z_face_fraction(i,j,k)*limited(p)%z(c,i,j,k)))/volume
        enddo
        if(present(internal_energy_rate)) then
          held_increment=( &
            areas(1)*(g%x_face_fraction(i-1,j,k)*held_delta%x(1,i-1,j,k)- &
              g%x_face_fraction(i,j,k)*held_delta%x(1,i,j,k))+ &
            areas(2)*(g%y_face_fraction(i,j-1,k)*held_delta%y(1,i,j-1,k)- &
              g%y_face_fraction(i,j,k)*held_delta%y(1,i,j,k))+ &
            areas(3)*(g%z_face_fraction(i,j,k-1)*held_delta%z(1,i,j,k-1)- &
              g%z_face_fraction(i,j,k)*held_delta%z(1,i,j,k)))/volume
          if(held_increment/=0) held_candidate%fields(p)%values(1,i,j,k)= &
            internal_energy_rate%fields(p)%values(1,i,j,k)+held_increment
          if(.not.ieee_is_finite(held_candidate%fields(p)%values(1,i,j,k))) then
            return
          endif
        endif
        endpoint=tree%nodes(p)%state(:,i,j,k)+dt*(candidate%fields(p)%values(:,i,j,k)+source%fields(p)%values(:,i,j,k))
        if(any(.not.ieee_is_finite(endpoint)).or.endpoint(1)<=0.or.any(endpoint(6:)<0)) return
      enddo
      enddo
      enddo
      end associate
    enddo
    do p=1,size(tree%nodes)
      if(depth>=0.and.tree%topology%nodes(p)%level/=depth) cycle
      call move_alloc(limited(p)%x,flux(p)%x);call move_alloc(limited(p)%y,flux(p)%y)
      call move_alloc(limited(p)%z,flux(p)%z)
      if(present(internal_energy_rate)) then
        call move_alloc(held_limited(p)%x,internal_energy_flux(p)%x)
        call move_alloc(held_limited(p)%y,internal_energy_flux(p)%y)
        call move_alloc(held_limited(p)%z,internal_energy_flux(p)%z)
      endif
    enddo
    if(present(internal_energy_rate)) call move_alloc(held_candidate%fields,internal_energy_rate%fields)
    call move_alloc(candidate%fields,hydro_rate%fields);ok=.true.
  end subroutine apply_eb_tree_flux_correction_3d

  subroutine compute_godunov_state_redist_eb_tree_hydro_3d(species,tree,predictor,extended_geometry,dt,hydro_rate,flux,ok, &
      difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
      target_depth,ancestor_start,ancestor_end,ancestor_fraction,source_start,source_end,source_fraction,baseline_flux, &
      spatial_rate_halos,old_state_halos, &
      registers,max_order,target_volume_fraction,correct_forward_flux, &
      internal_energy_state,internal_energy_rate,internal_energy_flux,internal_energy_registers, &
      internal_energy_baseline_flux,internal_energy_spatial_rate_halos,internal_energy_old_state_halos, &
      ancestor_internal_energy_start,ancestor_internal_energy_end,internal_energy_pressure_work,spatial_flux_weight)
    real(dp), optional, intent(in) :: spatial_flux_weight
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy_state,ancestor_internal_energy_start, &
      ancestor_internal_energy_end
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy_rate,internal_energy_spatial_rate_halos, &
      internal_energy_old_state_halos,internal_energy_pressure_work
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: internal_energy_flux(:),internal_energy_baseline_flux(:)
    type(state_redist_register_3d), allocatable, optional, intent(inout) :: internal_energy_registers(:)
    type(amr_tree_state_3d) :: pending_held_rate,pending_held_source,pending_held_old,pending_held_work
    type(eb_integrated_flux_3d), allocatable :: pending_held_flux(:),pending_held_baseline(:)
    type(state_redist_register_3d), allocatable :: pending_held_registers(:)
    integer :: held_options
    logical :: held,held_valid
    type(amr_tree_state_3d), optional, intent(inout) :: spatial_rate_halos,old_state_halos
    logical, optional, intent(in) :: correct_forward_flux
    type(nasa7_species), intent(in) :: species(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: predictor
    type(eb_geometry_3d), intent(in) :: extended_geometry(:)
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: hydro_rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: baseline_flux(:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: target_depth
    type(eb_reactive_tree_3d), optional, intent(in) :: ancestor_start,ancestor_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    type(amr_tree_state_3d), optional, intent(in) :: source_start,source_end
    real(dp), optional, intent(in) :: source_fraction(:)
    integer :: selected_depth
    real(dp), optional, intent(in) :: difmag,inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(state_redist_register_3d), allocatable, intent(inout) :: registers(:)
    integer, optional, intent(in) :: max_order
    real(dp), optional, intent(in) :: target_volume_fraction
    type(amr_tree_state_3d) :: pending_rate,pending_source,pending_old
    type(eb_integrated_flux_3d), allocatable :: pending_flux(:),pending_baseline(:)
    type(state_redist_register_3d), allocatable :: pending_registers(:)
    type(eb_geometry_3d), allocatable :: redist_geometry(:)
    integer :: p,n(3),nv,first_depth,last_depth,depth,padding
    logical :: good
    ok=.false.;selected_depth=-1
    if(present(target_depth)) selected_depth=target_depth
    held_options=count([present(internal_energy_state),present(internal_energy_rate), &
      present(internal_energy_flux),present(internal_energy_registers)])
    held=held_options==4
    held_valid=held_options==0.or.held_options==4
    held_valid=held_valid.and.(held.or..not.(present(internal_energy_baseline_flux).or. &
      present(internal_energy_spatial_rate_halos).or.present(internal_energy_old_state_halos).or. &
      present(internal_energy_pressure_work).or.present(ancestor_internal_energy_start).or. &
      present(ancestor_internal_energy_end)))
    held_valid=held_valid.and.(present(ancestor_internal_energy_start).eqv.present(ancestor_internal_energy_end))
    if(held) held_valid=held_valid.and.(present(ancestor_fraction).eqv.present(ancestor_internal_energy_start))
    if(.not.held_valid) return
    if(present(baseline_flux)) then
      if(held) then
      call compute_sdc_eb_tree_hydro_3d(species=species,tree=tree,predictor=predictor, &
        extended_geometry=extended_geometry,dt=dt,hydro_rate=pending_rate,flux=pending_flux,ok=good, &
        difmag=difmag,boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
        inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,target_depth=target_depth, &
        ancestor_start=ancestor_start,ancestor_end=ancestor_end,ancestor_fraction=ancestor_fraction, &
        source_start=source_start,source_end=source_end,source_fraction=source_fraction, &
        baseline_flux=pending_baseline,spatial_rate_halos=pending_source,old_state_halos=pending_old, &
        correct_forward_flux=correct_forward_flux, &
        internal_energy_state=internal_energy_state,internal_energy_rate=pending_held_rate, &
        internal_energy_flux=pending_held_flux,internal_energy_baseline_flux=pending_held_baseline, &
        internal_energy_spatial_rate_halos=pending_held_source,internal_energy_old_state_halos=pending_held_old, &
        ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end,internal_energy_pressure_work=pending_held_work, &
        spatial_flux_weight=spatial_flux_weight)
      else
      call compute_sdc_eb_tree_hydro_3d(species=species,tree=tree,predictor=predictor, &
        extended_geometry=extended_geometry,dt=dt,hydro_rate=pending_rate,flux=pending_flux,ok=good, &
        difmag=difmag,boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
        inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,target_depth=target_depth, &
        ancestor_start=ancestor_start,ancestor_end=ancestor_end,ancestor_fraction=ancestor_fraction, &
        source_start=source_start,source_end=source_end,source_fraction=source_fraction, &
        baseline_flux=pending_baseline,spatial_rate_halos=pending_source,old_state_halos=pending_old, &
        correct_forward_flux=correct_forward_flux, &
        spatial_flux_weight=spatial_flux_weight)
      endif
    else
      if(held) then
      call compute_sdc_eb_tree_hydro_3d(species=species,tree=tree,predictor=predictor, &
        extended_geometry=extended_geometry,dt=dt,hydro_rate=pending_rate,flux=pending_flux,ok=good, &
        difmag=difmag,boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
        inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,target_depth=target_depth, &
        ancestor_start=ancestor_start,ancestor_end=ancestor_end,ancestor_fraction=ancestor_fraction, &
        source_start=source_start,source_end=source_end,source_fraction=source_fraction, &
        spatial_rate_halos=pending_source,old_state_halos=pending_old, &
        correct_forward_flux=correct_forward_flux, &
        internal_energy_state=internal_energy_state,internal_energy_rate=pending_held_rate, &
        internal_energy_flux=pending_held_flux,internal_energy_baseline_flux=pending_held_baseline, &
        internal_energy_spatial_rate_halos=pending_held_source,internal_energy_old_state_halos=pending_held_old, &
        ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end,internal_energy_pressure_work=pending_held_work, &
        spatial_flux_weight=spatial_flux_weight)
      else
      call compute_sdc_eb_tree_hydro_3d(species=species,tree=tree,predictor=predictor, &
        extended_geometry=extended_geometry,dt=dt,hydro_rate=pending_rate,flux=pending_flux,ok=good, &
        difmag=difmag,boundary_faces=boundary_faces,inflow_state=inflow_state,inflow_mach=inflow_mach, &
        inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure,target_depth=target_depth, &
        ancestor_start=ancestor_start,ancestor_end=ancestor_end,ancestor_fraction=ancestor_fraction, &
        source_start=source_start,source_end=source_end,source_fraction=source_fraction, &
        spatial_rate_halos=pending_source,old_state_halos=pending_old, &
        correct_forward_flux=correct_forward_flux, &
        spatial_flux_weight=spatial_flux_weight)
      endif
    endif
    if(.not.good) return
    nv=reactive_nvar(size(species));first_depth=0
    if(present(target_depth)) first_depth=target_depth
    last_depth=first_depth
    if(.not.present(target_depth)) last_depth=tree%topology%level_count()-1
    allocate(redist_geometry(size(tree%nodes)))
    do p=1,size(tree%nodes)
      n=tree%topology%nodes(p)%shape
      padding=(extended_geometry(p)%nx-n(1))/2
      call extract_eb_geometry_3d(extended_geometry(p),[padding-5,padding-5,padding-5], &
        n+padding+6,redist_geometry(p),good)
      if(.not.good) exit
    enddo
    if(.not.good) return
    ! Pinned Hydro.cpp uses factor1 for ML StateRedist. The caller owns
    ! provisional coarse accumulation and final-iteration fine/flux acceptance.
    do depth=first_depth,last_depth
      call compute_state_redist_level_stage_3d(tree%topology,redist_geometry,pending_old,pending_source, &
        nv,depth,dt,pending_rate,pending_registers,good,boundary_faces,target_volume_fraction, &
        max_order,register_scale=1._dp)
      if(.not.good) return
      if(held) then
        call compute_state_redist_level_stage_3d(tree%topology,redist_geometry, &
          pending_held_old,pending_held_source,1,depth,dt,pending_held_rate,pending_held_registers,good, &
          boundary_faces,target_volume_fraction,max_order,register_scale=1._dp)
        if(.not.good) return
      endif
    enddo
    do p=1,size(tree%nodes)
      if(selected_depth>=0.and.tree%topology%nodes(p)%level/=selected_depth) cycle
      call post_eb_spatial_density_rate_3d(tree%nodes(p)%geometry,pending_rate%fields(p)%values,good)
      if(.not.good) return
    enddo
    if(held) then
      call move_alloc(pending_held_rate%fields,internal_energy_rate%fields)
      call move_alloc(pending_held_flux,internal_energy_flux)
      call move_alloc(pending_held_registers,internal_energy_registers)
      if(present(internal_energy_baseline_flux)) call move_alloc(pending_held_baseline,internal_energy_baseline_flux)
      if(present(internal_energy_spatial_rate_halos)) &
        call move_alloc(pending_held_source%fields,internal_energy_spatial_rate_halos%fields)
      if(present(internal_energy_old_state_halos)) &
        call move_alloc(pending_held_old%fields,internal_energy_old_state_halos%fields)
      if(present(internal_energy_pressure_work)) call move_alloc(pending_held_work%fields,internal_energy_pressure_work%fields)
    endif
    call move_alloc(pending_rate%fields,hydro_rate%fields)
    call move_alloc(pending_flux,flux)
    call move_alloc(pending_registers,registers)
    if(present(baseline_flux)) call move_alloc(pending_baseline,baseline_flux)
    ! These remain the raw spatial RHS and actual old state used by StateRedist.
    if(present(spatial_rate_halos)) call move_alloc(pending_source%fields,spatial_rate_halos%fields)
    if(present(old_state_halos)) call move_alloc(pending_old%fields,old_state_halos%fields)
    ok=.true.
  end subroutine compute_godunov_state_redist_eb_tree_hydro_3d

  subroutine compute_sdc_eb_tree_hydro_3d(species,tree,predictor,extended_geometry,dt,hydro_rate,flux,ok, &
      difmag,boundary_faces,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
      target_depth,ancestor_start,ancestor_end,ancestor_fraction,source_start,source_end,source_fraction,baseline_flux, &
      spatial_rate_halos,old_state_halos,correct_forward_flux, &
      internal_energy_state,internal_energy_rate,internal_energy_flux,internal_energy_baseline_flux, &
      internal_energy_spatial_rate_halos,internal_energy_old_state_halos, &
      ancestor_internal_energy_start,ancestor_internal_energy_end,internal_energy_pressure_work,spatial_flux_weight)
    real(dp), optional, intent(in) :: spatial_flux_weight
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy_state,ancestor_internal_energy_start, &
      ancestor_internal_energy_end
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy_rate,internal_energy_spatial_rate_halos, &
      internal_energy_old_state_halos,internal_energy_pressure_work
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: internal_energy_flux(:),internal_energy_baseline_flux(:)
    type(amr_tree_state_3d) :: held_rates,held_spatial,held_old,held_work
    type(eb_integrated_flux_3d), allocatable :: held_faces(:),held_baselines(:),held_local_baseline
    type(eb_integrated_flux_3d) :: held_full_flux
    real(dp), allocatable :: held_full_rate(:,:,:),held_full_work(:,:,:)
    integer :: held_options
    logical :: held,held_valid
    type(amr_tree_state_3d), optional, intent(inout) :: spatial_rate_halos,old_state_halos
    logical, optional, intent(in) :: correct_forward_flux
    type(nasa7_species), intent(in) :: species(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: predictor
    type(eb_geometry_3d), intent(in) :: extended_geometry(:)
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: hydro_rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    type(eb_integrated_flux_3d), allocatable, optional, intent(inout) :: baseline_flux(:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: target_depth
    type(eb_reactive_tree_3d), optional, intent(in) :: ancestor_start,ancestor_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    type(amr_tree_state_3d), optional, intent(in) :: source_start,source_end
    real(dp), optional, intent(in) :: source_fraction(:)
    integer :: selected_depth
    real(dp), optional, intent(in) :: difmag,inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(amr_tree_state_3d) :: rates,pending_halos,pending_old
    type(eb_integrated_flux_3d), allocatable :: candidate_flux(:)
    type(eb_integrated_flux_3d) :: full_flux
    real(dp), allocatable :: held_ghost(:,:,:)
    type(eb_integrated_flux_3d), allocatable :: candidate_baseline(:),local_baseline
    type(eb_geometry_3d), allocatable :: halo_geometry(:)
    type(eb_geometry_3d) :: working,center,correction_geometry
    type(eb_plm_boundary_3d) :: physical_boundary
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),vf(:,:,:),flat(:,:,:),source(:,:,:,:),full_rate(:,:,:,:)
    real(dp) :: coefficient,boundary(reactive_nvar(size(species))),wall(reactive_nvar(size(species)))
    character(len=32) :: faces(6)
    integer :: nv,node,n(3),level,width,grow,available_padding,offset_geometry
    logical :: valid,corrected
    ok=.false.;nv=reactive_nvar(size(species))
    held_options=count([present(internal_energy_state),present(internal_energy_rate),present(internal_energy_flux)])
    held=held_options==3
    held_valid=held_options==0.or.held_options==3
    held_valid=held_valid.and.(held.or..not.(present(internal_energy_baseline_flux).or. &
      present(internal_energy_spatial_rate_halos).or.present(internal_energy_old_state_halos).or. &
      present(internal_energy_pressure_work).or.present(ancestor_internal_energy_start).or. &
      present(ancestor_internal_energy_end)))
    held_valid=held_valid.and.(present(ancestor_internal_energy_start).eqv.present(ancestor_internal_energy_end))
    if(held) held_valid=held_valid.and.(present(ancestor_fraction).eqv.present(ancestor_internal_energy_start))
    if(.not.held_valid) return
    corrected=.false.
    if(present(correct_forward_flux)) corrected=correct_forward_flux
    if(.not.ieee_is_finite(dt).or.dt<=0) return
    if(.not.valid_tree(tree,nv)) return
    if(held) then
      valid=stored_field_valid(internal_energy_state)
      if(present(ancestor_internal_energy_start)) valid=valid.and. &
        stored_field_valid(ancestor_internal_energy_start).and.stored_field_valid(ancestor_internal_energy_end)
      if(.not.valid) return
    endif
    selected_depth=-1
    if(present(target_depth)) then
      selected_depth=target_depth
      if(selected_depth<0.or.selected_depth>=tree%topology%level_count()) return
    end if
    if(present(ancestor_start).or.present(ancestor_end).or.present(ancestor_fraction)) then
      if(selected_depth<0) return
      if(.not.(present(ancestor_start).and.present(ancestor_end).and.present(ancestor_fraction))) return
    end if
    if(present(source_start).or.present(source_end).or.present(source_fraction)) then
      if(selected_depth<0) return
      if(.not.(present(source_start).and.present(source_end).and.present(source_fraction))) return
    end if

    width=6
    if(present(spatial_rate_halos).or.present(old_state_halos).or.present(internal_energy_spatial_rate_halos).or. &
      present(internal_energy_old_state_halos)) width=9
    if(corrected) width=10
    grow=width-6
    if(size(extended_geometry)/=size(tree%nodes)) return
    allocate(halo_geometry(size(tree%nodes)))
    do node=1,size(tree%nodes)
      valid=extended_geometry(node)%is_valid()
      if(.not.valid) return
      n=tree%topology%nodes(node)%shape;available_padding=(extended_geometry(node)%nx-n(1))/2
      valid=available_padding>=width.and.all([extended_geometry(node)%nx, &
        extended_geometry(node)%ny,extended_geometry(node)%nz]==n+2*available_padding)
      if(.not.valid) return
      offset_geometry=available_padding-width
      call extract_eb_geometry_3d(extended_geometry(node),[1,1,1]+offset_geometry, &
        n+2*width+offset_geometry,halo_geometry(node),valid)
      if(.not.valid) return
    enddo
    if(.not.allocated(predictor%fields)) return
    if(size(predictor%fields)/=size(tree%nodes)) return
    coefficient=.1_dp
    if(present(difmag)) coefficient=difmag
    if(.not.ieee_is_finite(coefficient).or.coefficient<0) return
    faces='outflow'
    if(present(boundary_faces)) faces=boundary_faces
    if(present(inflow_mach)) then
      where(inflow_mach>=0) faces='inflow'
    end if
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if(.not.valid) return
    if(any(domain>int(huge(1)-4,int64))) return
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      if(.not.allocated(predictor%fields(node)%values)) return
      if(any(shape(predictor%fields(node)%values)/=[nv,n])) return
      if(any([halo_geometry(node)%nx,halo_geometry(node)%ny,halo_geometry(node)%nz]/=n+2*width)) return
      call extract_eb_geometry_3d(halo_geometry(node),[1,1,1]+width,n+width,center,valid)
      if(.not.valid) return
      if(.not.matching_eb_metrics_3d(center,tree%nodes(node)%geometry)) return
    end do
    allocate(rates%fields(size(tree%nodes)),candidate_flux(size(tree%nodes)))
    if(width>=9) allocate(pending_halos%fields(size(tree%nodes)),pending_old%fields(size(tree%nodes)))
    if(present(baseline_flux)) allocate(candidate_baseline(size(tree%nodes)))
    if(present(baseline_flux).or.corrected.or.present(internal_energy_baseline_flux)) allocate(local_baseline)
    if(held) then
      allocate(held_rates%fields(size(tree%nodes)),held_faces(size(tree%nodes)))

      if(present(internal_energy_baseline_flux)) then
        allocate(held_baselines(size(tree%nodes)))

      endif
      if(present(internal_energy_baseline_flux).or.corrected) then
        allocate(held_local_baseline )

      endif
      if(width>=9) then
        allocate(held_spatial%fields(size(tree%nodes)),held_old%fields(size(tree%nodes)))

      endif
      if(present(internal_energy_pressure_work)) then
        allocate(held_work%fields(size(tree%nodes)))

      endif
    endif
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      allocate(rates%fields(node)%values(nv,n(1),n(2),n(3)),source=0._dp)
      if(held) then
        allocate(held_rates%fields(node)%values(1,n(1),n(2),n(3)),source=0._dp)
        if(present(internal_energy_pressure_work)) &
          allocate(held_work%fields(node)%values(1,n(1),n(2),n(3)),source=0._dp)
      endif
      if(selected_depth>=0.and.tree%topology%nodes(node)%level/=selected_depth) cycle
      allocate(ghost(nv,1-width:n(1)+width,1-width:n(2)+width,1-width:n(3)+width), &
        gt(1-width:n(1)+width,1-width:n(2)+width,1-width:n(3)+width), &
        vf(1-width:n(1)+width,1-width:n(2)+width,1-width:n(3)+width), &
        flat(0:n(1)+2*width-5,0:n(2)+2*width-5,0:n(3)+2*width-5), &
        source(nv,n(1)+2*width-6,n(2)+2*width-6,n(3)+2*width-6), &
        full_rate(nv,n(1)+2*grow,n(2)+2*grow,n(3)+2*grow))
      if(held) then
      allocate(held_ghost(1-width:n(1)+width,1-width:n(2)+width,1-width:n(3)+width))
      call fill_reactive_eb_tree_halo_3d(species,tree,node,ghost,gt,valid, &
        inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,ghost_volume_fraction=vf,halo_width=width, &
        target_volume_fraction=halo_geometry(node)%volume_fraction,restrict_snapshot=.false., &
        ancestor_start=ancestor_start,ancestor_end=ancestor_end,ancestor_fraction=ancestor_fraction, &
        internal_energy=internal_energy_state,ghost_internal_energy=held_ghost, &
        ancestor_internal_energy_start=ancestor_internal_energy_start, &
        ancestor_internal_energy_end=ancestor_internal_energy_end)
      else
      call fill_reactive_eb_tree_halo_3d(species,tree,node,ghost,gt,valid, &
        inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,ghost_volume_fraction=vf,halo_width=width, &
        target_volume_fraction=halo_geometry(node)%volume_fraction,restrict_snapshot=.false., &
        ancestor_start=ancestor_start,ancestor_end=ancestor_end,ancestor_fraction=ancestor_fraction)
      endif
      if(.not.valid) return
      call extract_eb_geometry_3d(halo_geometry(node),[4,4,4],n+2*width-3,working,valid)
      if(.not.valid) return
      call fill_sdc_eb_tree_source_halo_3d(tree,predictor,node,width-3,source,valid, &
        boundary_faces=faces,target_volume_fraction=working%volume_fraction, &
        ancestor_start=source_start,ancestor_end=source_end,ancestor_fraction=source_fraction)
      if(.not.valid) return
      call compute_reactive_eb_flattening_3d(species,ghost,gt,vf,n+2*width-6,flat,valid)
      if(.not.valid) return
      level=tree%topology%nodes(node)%level+1
      physical_boundary%lower=width-2-int(origin(:,node))
      physical_boundary%upper=width-3+int(domain(:,level)-origin(:,node))
      physical_boundary%faces=faces
      if(held) then
        allocate(held_full_rate(n(1)+2*grow,n(2)+2*grow,n(3)+2*grow), &
          held_full_work(n(1)+2*grow,n(2)+2*grow,n(3)+2*grow))
      call compute_sdc_eb_hydro_3d(species,working,ghost,gt,vf,dt,flat,source,full_rate, &
        full_flux,boundary,wall,valid,difmag=coefficient,boundary=physical_boundary,baseline_flux=local_baseline, &
        internal_energy_state=held_ghost,internal_energy_rate=held_full_rate, &
        internal_energy_integrated_flux=held_full_flux,internal_energy_baseline_flux=held_local_baseline, &
        internal_energy_pressure_work=held_full_work, &
        spatial_flux_weight=spatial_flux_weight)
      else
      call compute_sdc_eb_hydro_3d(species,working,ghost,gt,vf,dt,flat,source,full_rate, &
        full_flux,boundary,wall,valid,difmag=coefficient,boundary=physical_boundary,baseline_flux=local_baseline, &
        spatial_flux_weight=spatial_flux_weight)
      endif
      if(.not.valid) return
      if(corrected) then
        call extract_eb_geometry_3d(working,[4,4,4],n+2*grow+3,correction_geometry,valid)
        if(.not.valid) return
        if(held) then
        call correct_reactive_eb_spatial_hydro_flux_3d(correction_geometry, &
          ghost(:,1-grow:n(1)+grow,1-grow:n(2)+grow,1-grow:n(3)+grow), &
          source(:,4:n(1)+2*grow+3,4:n(2)+2*grow+3,4:n(3)+2*grow+3), &
          dt,full_rate,full_flux,local_baseline,valid, &
          internal_energy_rate=held_full_rate,internal_energy_flux=held_full_flux, &
          internal_energy_baseline=held_local_baseline)
        else
        call correct_reactive_eb_spatial_hydro_flux_3d(correction_geometry, &
          ghost(:,1-grow:n(1)+grow,1-grow:n(2)+grow,1-grow:n(3)+grow), &
          source(:,4:n(1)+2*grow+3,4:n(2)+2*grow+3,4:n(3)+2*grow+3), &
          dt,full_rate,full_flux,local_baseline,valid)
        endif
        if(.not.valid) return
      endif
      rates%fields(node)%values=full_rate(:,grow+1:n(1)+grow,grow+1:n(2)+grow,grow+1:n(3)+grow)
      call crop_flux(full_flux,candidate_flux(node),n,grow)
      if(present(baseline_flux)) call crop_flux(local_baseline,candidate_baseline(node),n,grow)
      if(held) then
        held_rates%fields(node)%values(1,:,:,:)= &
          held_full_rate(grow+1:n(1)+grow,grow+1:n(2)+grow,grow+1:n(3)+grow)
        call crop_flux(held_full_flux,held_faces(node),n,grow)

        if(present(internal_energy_baseline_flux)) then
          call crop_flux(held_local_baseline,held_baselines(node),n,grow)

        endif
        if(present(internal_energy_pressure_work)) held_work%fields(node)%values(1,:,:,:)= &
          held_full_work(grow+1:n(1)+grow,grow+1:n(2)+grow,grow+1:n(3)+grow)
        if(width>=9) then
          allocate(held_spatial%fields(node)%values(1,n(1)+12,n(2)+12,n(3)+12), &
            held_old%fields(node)%values(1,n(1)+12,n(2)+12,n(3)+12))

          held_spatial%fields(node)%values=0
          held_spatial%fields(node)%values(1,4:n(1)+9,4:n(2)+9,4:n(3)+9)= &
            held_full_rate(grow-2:n(1)+grow+3,grow-2:n(2)+grow+3,grow-2:n(3)+grow+3)
          held_old%fields(node)%values(1,:,:,:)=held_ghost(1-width+grow:n(1)+width-grow, &
            1-width+grow:n(2)+width-grow,1-width+grow:n(3)+width-grow)
        endif
        deallocate(held_full_rate,held_full_work,held_ghost)
      endif
      if(width>=9) then
        allocate(pending_halos%fields(node)%values(nv,n(1)+12,n(2)+12,n(3)+12),source=0._dp)
        pending_halos%fields(node)%values(:,4:n(1)+9,4:n(2)+9,4:n(3)+9)= &
          full_rate(:,grow-2:n(1)+grow+3,grow-2:n(2)+grow+3,grow-2:n(3)+grow+3)
        pending_old%fields(node)%values=ghost(:,1-width+grow:n(1)+width-grow,1-width+grow:n(2)+width-grow, &
          1-width+grow:n(3)+width-grow)
      endif
      deallocate(ghost,gt,vf,flat,source,full_rate)
    end do
    if(width>=9) then
      call refresh_peer_source_halos(valid)
      if(.not.valid) return
    endif
    if(present(baseline_flux)) call move_alloc(candidate_baseline,baseline_flux)
    if(present(spatial_rate_halos)) call move_alloc(pending_halos%fields,spatial_rate_halos%fields)
    if(present(old_state_halos)) call move_alloc(pending_old%fields,old_state_halos%fields)
    if(held) then
      call move_alloc(held_rates%fields,internal_energy_rate%fields)
      call move_alloc(held_faces,internal_energy_flux)
      if(present(internal_energy_baseline_flux)) call move_alloc(held_baselines,internal_energy_baseline_flux)
      if(present(internal_energy_spatial_rate_halos)) &
        call move_alloc(held_spatial%fields,internal_energy_spatial_rate_halos%fields)
      if(present(internal_energy_old_state_halos)) call move_alloc(held_old%fields,internal_energy_old_state_halos%fields)
      if(present(internal_energy_pressure_work)) call move_alloc(held_work%fields,internal_energy_pressure_work%fields)
    endif
    hydro_rate=rates;call move_alloc(candidate_flux,flux);ok=.true.
  contains
    logical function stored_field_valid(field) result(good)
      type(amr_tree_state_3d), intent(in) :: field
      integer :: p,d(3)
      good=.false.
      if(.not.allocated(field%fields)) return
      if(size(field%fields)/=size(tree%nodes).or.lbound(field%fields,1)/=1) return
      do p=1,size(tree%nodes)

        if(.not.allocated(field%fields(p)%values)) return
        d=tree%topology%nodes(p)%shape
        if(any(shape(field%fields(p)%values)/=[1,d]).or.any(lbound(field%fields(p)%values)/=1)) return
        if(any(.not.ieee_is_finite(field%fields(p)%values(1,:,:,:)).and.tree%nodes(p)%geometry%volume_fraction>0)) return
      enddo
      good=.true.
    end function stored_field_valid
    subroutine crop_flux(full,cropped,n,g)
      type(eb_integrated_flux_3d), intent(in) :: full
      type(eb_integrated_flux_3d), intent(out) :: cropped
      integer, intent(in) :: n(3),g
      integer :: nc
      nc=size(full%x,1)
      allocate(cropped%x(nc,0:n(1),n(2),n(3)),cropped%y(nc,n(1),0:n(2),n(3)), &
        cropped%z(nc,n(1),n(2),0:n(3)))
      cropped%x=full%x(:,g:n(1)+g,g+1:n(2)+g,g+1:n(3)+g)
      cropped%y=full%y(:,g+1:n(1)+g,g:n(2)+g,g+1:n(3)+g)
      cropped%z=full%z(:,g+1:n(1)+g,g+1:n(2)+g,g:n(3)+g)
    end subroutine crop_flux
    subroutine refresh_peer_source_halos(success)
      logical, intent(out) :: success
      integer(int64), allocatable :: origin(:,:),domain(:,:)
      integer, allocatable :: ratios(:)
      integer(int64) :: cell(3)
      integer :: p,q,level,i,j,k,index(3),a,m(3)
      character(len=32) :: faces(6)
      logical :: inside
      faces='outflow'
      if(present(boundary_faces)) faces=boundary_faces
      call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,success)
      if(.not.success) return
      do p=1,size(tree%nodes)
        if(.not.allocated(pending_halos%fields(p)%values)) cycle
        level=tree%topology%nodes(p)%level;m=tree%topology%nodes(p)%shape
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
          do q=1,size(tree%nodes)
            if(tree%topology%nodes(q)%level/=level) cycle
            if(any(cell<origin(:,q)).or.any(cell>=origin(:,q)+int(tree%topology%nodes(q)%shape,int64))) cycle
            index=int(cell-origin(:,q))+1
            pending_halos%fields(p)%values(:,i,j,k)=rates%fields(q)%values(:,index(1),index(2),index(3))
            if(held) held_spatial%fields(p)%values(:,i,j,k)=held_rates%fields(q)%values(:,index(1),index(2),index(3))
            exit
          enddo
        enddo
        enddo
        enddo
      enddo
    end subroutine refresh_peer_source_halos
  end subroutine compute_sdc_eb_tree_hydro_3d

  subroutine advance_reactive_eb_tree_hydro_3d(species,solver,dt,tree,boundary_change,wall_change,ok, &
      inflow_mach,inflow_temperature,inflow_pressure,difmag,extended_geometry,boundary_faces,inflow_state)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry(:)
    type(eb_geometry_3d) :: working,center
    type(eb_plm_boundary_3d) :: physical_boundary
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: flat(:,:,:)
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    type(eb_reactive_tree_3d), intent(inout) :: tree
    real(dp), intent(inout) :: boundary_change(:),wall_change(:)
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: difmag
    real(dp) :: coefficient
    real(dp), allocatable :: vf(:,:,:)
    logical, intent(out) :: ok
    type(eb_reactive_tree_3d) :: snapshot,candidate
    type(eb_integrated_flux_3d), allocatable :: flux(:)
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),source(:,:,:,:),boundary(:),wall(:)
    logical, allocatable :: exposed(:,:,:)
    integer :: nv,node,child,n(3),lo(3),hi(3),c,level
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if (.not.ieee_is_finite(dt).or.dt<=0._dp) return
    if (size(boundary_change)/=nv.or.size(wall_change)/=nv) return
    if (.not.valid_tree(tree,nv)) return
    coefficient=.1_dp
    if(present(difmag)) coefficient=difmag
    if(.not.ieee_is_finite(coefficient).or.coefficient<0) return
    if(present(extended_geometry)) then
      if(solver/='pelec'.or.size(extended_geometry)/=size(tree%nodes)) return
      call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
      if(.not.valid) return
      if(any(domain>int(huge(1)-4,int64))) return
      do node=1,size(tree%nodes)
        n=tree%topology%nodes(node)%shape
        if(any([extended_geometry(node)%nx,extended_geometry(node)%ny,extended_geometry(node)%nz]/=n+12)) return
        call extract_eb_geometry_3d(extended_geometry(node),[7,7,7],n+6,center,valid)
        if(.not.valid) return
        if(.not.matching_eb_metrics_3d(center,tree%nodes(node)%geometry)) return
      end do
    end if
    snapshot=tree
    call average_down_reactive_eb_tree_3d(species,snapshot,valid)
    if(.not.valid) return
    candidate=snapshot
    allocate(flux(size(tree%nodes)),boundary(nv),wall(nv));wall=0._dp
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      if(present(extended_geometry)) then
        allocate(ghost(nv,-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),gt(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), &
          vf(-5:n(1)+6,-5:n(2)+6,-5:n(3)+6),flat(0:n(1)+7,0:n(2)+7,0:n(3)+7))
        call fill_reactive_eb_tree_halo_3d(species,snapshot,node,ghost,gt,valid, &
          inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
          boundary_faces=boundary_faces,inflow_state=inflow_state, &
          ghost_volume_fraction=vf,halo_width=6,target_volume_fraction=extended_geometry(node)%volume_fraction)
        if(.not.valid) return
        call extract_eb_geometry_3d(extended_geometry(node),[4,4,4],n+9,working,valid)
        if(.not.valid) return
        call compute_reactive_eb_flattening_3d(species,ghost,gt,vf,n+6,flat,valid)
        if(.not.valid) return
        level=tree%topology%nodes(node)%level+1
        physical_boundary%lower=4-int(origin(:,node))
        physical_boundary%upper=3+int(domain(:,level)-origin(:,node))
        physical_boundary%faces='outflow'
        if(present(boundary_faces)) physical_boundary%faces=boundary_faces
        if(present(inflow_mach)) then
          where(inflow_mach>=0) physical_boundary%faces='inflow'
        end if
        call advance_reactive_eb_plm_ctu_3d(species,working,ghost,gt,vf,dt,flat,candidate%nodes(node)%state, &
          candidate%nodes(node)%temperature,valid,difmag=coefficient,integrated_flux=flux(node),boundary=physical_boundary)
        if(.not.valid) return
        do c=1,nv
          where(snapshot%nodes(node)%geometry%volume_fraction<=0) &
            candidate%nodes(node)%state(c,:,:,:)=snapshot%nodes(node)%state(c,:,:,:)
        end do
        where(snapshot%nodes(node)%geometry%volume_fraction<=0) &
          candidate%nodes(node)%temperature=snapshot%nodes(node)%temperature
        deallocate(flat)
      else
      allocate(ghost(nv,0:n(1)+1,0:n(2)+1,0:n(3)+1),gt(0:n(1)+1,0:n(2)+1,0:n(3)+1))
      if(coefficient>0) allocate(vf(0:n(1)+1,0:n(2)+1,0:n(3)+1))
      call fill_reactive_eb_tree_halo_3d(species,snapshot,node,ghost,gt,valid, &
        inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,ghost_volume_fraction=vf)
      if(.not.valid) return
      call advance_reactive_eb_euler_3d(species,snapshot%nodes(node)%state,snapshot%nodes(node)%temperature, &
        snapshot%nodes(node)%geometry,solver,dt,candidate%nodes(node)%state,candidate%nodes(node)%temperature, &
        valid,ghost,gt,flux(node),difmag=coefficient,ghost_volume_fraction=vf)
      end if
      if(.not.valid) return
      deallocate(ghost,gt)
      if(allocated(vf)) deallocate(vf)
      allocate(source(nv,n(1),n(2),n(3)),exposed(n(1),n(2),n(3)));exposed=.true.
      do child=node+1,size(tree%nodes)
        if (tree%topology%nodes(child)%parent/=node) cycle
        call bounds(tree%topology,child,lo,hi)
        exposed(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
      end do
      call reactive_eb_slip_wall_source_3d(species,snapshot%nodes(node)%state,snapshot%nodes(node)%temperature, &
        snapshot%nodes(node)%geometry,source,valid)
      if(.not.valid) return
      associate(g=>snapshot%nodes(node)%geometry)
      do c=1,nv
        wall(c)=wall(c)+dt*g%dx*g%dy*g%dz*sum(source(c,:,:,:)*g%volume_fraction,mask=exposed)
      end do
      end associate
      deallocate(source,exposed)
    end do
    call synchronize_reactive_eb_tree_flux_3d(species,candidate,flux,valid,boundary_faces=boundary_faces)
    if(.not.valid) return
    n=tree%topology%nodes(1)%shape
    associate(g=>tree%nodes(1)%geometry,f=>flux(1))
    do c=1,nv
      boundary(c)=g%dy*g%dz*sum(g%x_face_fraction(0,:,:)*f%x(c,0,:,:)- &
        g%x_face_fraction(n(1),:,:)*f%x(c,n(1),:,:))+ &
        g%dx*g%dz*sum(g%y_face_fraction(:,0,:)*f%y(c,:,0,:)- &
        g%y_face_fraction(:,n(2),:)*f%y(c,:,n(2),:))+ &
        g%dx*g%dy*sum(g%z_face_fraction(:,:,0)*f%z(c,:,:,0)- &
        g%z_face_fraction(:,:,n(3))*f%z(c,:,:,n(3)))
    end do
    end associate
    if (.not.all(ieee_is_finite(boundary)).or..not.all(ieee_is_finite(wall))) return
    tree=candidate;boundary_change=boundary;wall_change=wall;ok=.true.
  end subroutine advance_reactive_eb_tree_hydro_3d

  ! Each level advances together, with ratio fine substeps, time-interpolated
  ! ancestor donors, per-substep peer fluxes and interval-local reflux.
  subroutine advance_reactive_eb_tree_subcycle_hydro_3d(species,solver,dt,tree,boundary_change,wall_change,ok, &
      inflow_mach,inflow_temperature,inflow_pressure,difmag,extended_geometry,boundary_faces,inflow_state)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry(:)
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    type(eb_reactive_tree_3d), intent(inout) :: tree
    real(dp), intent(inout) :: boundary_change(:),wall_change(:)
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: difmag
    logical, intent(out) :: ok
    call advance_reactive_eb_tree_subcycle_operator_3d(species,solver,dt,tree,boundary_change,wall_change,ok, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
      boundary_faces=boundary_faces,inflow_state=inflow_state, &
            difmag=difmag,extended_geometry=extended_geometry)
  end subroutine advance_reactive_eb_tree_subcycle_hydro_3d

  subroutine advance_reactive_eb_tree_subcycle_transport_3d(species,transport,dt,viscosity,conduction, &
      diffusion,barodiffusion,tree,boundary_change,minimum_theta,ok, &
      inflow_mach,inflow_temperature,inflow_pressure,boundary_faces,inflow_state,wall_temperature)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: dt
    logical, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    type(eb_reactive_tree_3d), intent(inout) :: tree
    real(dp), intent(inout) :: boundary_change(:),minimum_theta
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: wall_temperature(6)
    logical, intent(out) :: ok
    real(dp) :: wall(size(boundary_change))
    call advance_reactive_eb_tree_subcycle_operator_3d(species,'',dt,tree,boundary_change,wall,ok, &
      transport,viscosity,conduction,diffusion,barodiffusion,minimum_theta, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
      boundary_faces=boundary_faces,inflow_state=inflow_state,wall_temperature=wall_temperature)
  end subroutine advance_reactive_eb_tree_subcycle_transport_3d

  subroutine advance_reactive_eb_tree_subcycle_operator_3d(species,solver,dt,tree,boundary_change,wall_change,ok, &
      transport,viscosity,conduction,diffusion,barodiffusion,minimum_theta, &
      inflow_mach,inflow_temperature,inflow_pressure,difmag,extended_geometry,boundary_faces,inflow_state,wall_temperature)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry(:)
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    type(eb_reactive_tree_3d), intent(inout) :: tree
    real(dp), intent(inout) :: boundary_change(:),wall_change(:)
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: wall_temperature(6)
    real(dp), optional, intent(in) :: difmag
    logical, intent(out) :: ok
    type(gas_transport_species), optional, intent(in) :: transport(:)
    logical, optional, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), optional, intent(inout) :: minimum_theta
    type(eb_reactive_tree_3d) :: candidate,starts,ends
    type(eb_integrated_flux_3d), allocatable :: total(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: clock_start(:),clock_dt(:),wall(:),boundary(:)
    real(dp) :: theta,coefficient
    type(eb_geometry_3d) :: center
    integer :: nv,levels,c,n(3),node
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if(present(wall_temperature)) then
      if(.not.present(boundary_faces)) return
      if(.not.all(ieee_is_finite(wall_temperature)).or.any(wall_temperature<=0)) return
    end if
    if (.not.ieee_is_finite(dt).or.dt<=0._dp) return
    if (size(boundary_change)/=nv.or.size(wall_change)/=nv) return
    if (.not.valid_tree(tree,nv)) return
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if (.not.valid) return
    if(present(extended_geometry)) then
      if(present(transport).or.solver/='pelec') return
      if(size(extended_geometry)/=size(tree%nodes)) return
      if(any(domain>int(huge(1)-4,int64))) return
      do node=1,size(tree%nodes)
        n=tree%topology%nodes(node)%shape
        if(any([extended_geometry(node)%nx,extended_geometry(node)%ny,extended_geometry(node)%nz]/=n+12)) return
        call extract_eb_geometry_3d(extended_geometry(node),[7,7,7],n+6,center,valid)
        if(.not.valid) return
        if(.not.matching_eb_metrics_3d(center,tree%nodes(node)%geometry)) return
      end do
    end if
    if (present(transport)) then
      if (.not.(present(viscosity).and.present(conduction).and.present(diffusion).and. &
        present(barodiffusion).and.present(minimum_theta))) return
      if (any(ratios>huge(1)/ratios)) return
    end if
    coefficient=.1_dp
    if(present(transport)) coefficient=0
    if(present(difmag)) coefficient=difmag
    if(.not.ieee_is_finite(coefficient).or.coefficient<0) return
    if(present(transport).and.coefficient>0) return
    theta=1._dp
    levels=size(ratios);candidate=tree
    call average_down_reactive_eb_tree_3d(species,candidate,valid)
    if (.not.valid) return
    starts=candidate;ends=candidate
    allocate(clock_start(levels),clock_dt(levels),wall(nv),boundary(nv));wall=0._dp
    call advance_level(0,0._dp,dt,total,valid)
    if (.not.valid) return
    n=tree%topology%nodes(1)%shape
    associate(g=>tree%nodes(1)%geometry,f=>total(1))
    do c=1,nv
      boundary(c)=g%dy*g%dz*sum(g%x_face_fraction(0,:,:)*f%x(c,0,:,:)- &
        g%x_face_fraction(n(1),:,:)*f%x(c,n(1),:,:))+ &
        g%dx*g%dz*sum(g%y_face_fraction(:,0,:)*f%y(c,:,0,:)- &
        g%y_face_fraction(:,n(2),:)*f%y(c,:,n(2),:))+ &
        g%dx*g%dy*sum(g%z_face_fraction(:,:,0)*f%z(c,:,:,0)- &
        g%z_face_fraction(:,:,n(3))*f%z(c,:,:,n(3)))
    end do
    end associate
    if (.not.all(ieee_is_finite(boundary)).or..not.all(ieee_is_finite(wall))) return
    tree=candidate;boundary_change=boundary;wall_change=wall
    if (present(minimum_theta)) minimum_theta=theta
    ok=.true.
  contains
    recursive subroutine advance_level(level,time,interval,flux,good)
      integer, intent(in) :: level
      real(dp), intent(in) :: time,interval
      type(eb_integrated_flux_3d), allocatable, intent(out) :: flux(:)
      logical, intent(out) :: good
      type(eb_integrated_flux_3d), allocatable :: fine_flux(:)
      type(eb_integrated_flux_3d) :: owned
      real(dp) :: period(3),length(3)
      integer :: seam,tangent
      logical :: overlap
      type(eb_reactive_tree_3d) :: snapshot
      type(eb_geometry_3d) :: working
      type(eb_plm_boundary_3d) :: physical_boundary
      real(dp), allocatable :: flat(:,:,:)
      real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),source(:,:,:,:),fraction(:),vf(:,:,:)
      logical, allocatable :: exposed(:,:,:)
      real(dp) :: elapsed,sub_dt,stage_theta
      integer :: node,peer,owner,receiver,axis,child,substep,d(3),lo(3),hi(3),component,substeps
      integer(int64) :: offset(3)
      logical :: valid,faces(6)
      character(len=16) :: patch_faces(6)
      good=.false.;snapshot=candidate
      clock_start(level+1)=time;clock_dt(level+1)=interval
      allocate(flux(size(tree%nodes)),fraction(level))
      do node=1,size(tree%nodes)
        d=tree%topology%nodes(node)%shape
        allocate(flux(node)%x(nv,0:d(1),d(2),d(3)),flux(node)%y(nv,d(1),0:d(2),d(3)), &
          flux(node)%z(nv,d(1),d(2),0:d(3)))
        flux(node)%x=0._dp;flux(node)%y=0._dp;flux(node)%z=0._dp
        if (tree%topology%nodes(node)%level/=level) cycle
        starts%nodes(node)=snapshot%nodes(node)
        fraction=max(0._dp,min(1._dp,(time-clock_start(:level))/clock_dt(:level)))
        if(present(extended_geometry)) then
          allocate(ghost(nv,-5:d(1)+6,-5:d(2)+6,-5:d(3)+6),gt(-5:d(1)+6,-5:d(2)+6,-5:d(3)+6), &
            vf(-5:d(1)+6,-5:d(2)+6,-5:d(3)+6),flat(0:d(1)+7,0:d(2)+7,0:d(3)+7))
          call fill_reactive_eb_tree_halo_3d(species,snapshot,node,ghost,gt,valid,starts,ends,fraction, &
            inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
            boundary_faces=boundary_faces,inflow_state=inflow_state, &
            ghost_volume_fraction=vf,halo_width=6,target_volume_fraction=extended_geometry(node)%volume_fraction)
          if(.not.valid) return
          call extract_eb_geometry_3d(extended_geometry(node),[4,4,4],d+9,working,valid)
          if(.not.valid) return
          call compute_reactive_eb_flattening_3d(species,ghost,gt,vf,d+6,flat,valid)
          if(.not.valid) return
          physical_boundary%lower=4-int(origin(:,node))
          physical_boundary%upper=3+int(domain(:,level+1)-origin(:,node))
          physical_boundary%faces='outflow'
          if(present(boundary_faces)) physical_boundary%faces=boundary_faces
          if(present(inflow_mach)) then
            where(inflow_mach>=0) physical_boundary%faces='inflow'
          end if
        else
        allocate(ghost(nv,0:d(1)+1,0:d(2)+1,0:d(3)+1),gt(0:d(1)+1,0:d(2)+1,0:d(3)+1))
        if(coefficient>0) allocate(vf(0:d(1)+1,0:d(2)+1,0:d(3)+1))
        call fill_reactive_eb_tree_halo_3d(species,snapshot,node,ghost,gt,valid,starts,ends,fraction, &
          inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
          boundary_faces=boundary_faces,inflow_state=inflow_state,ghost_volume_fraction=vf)
        if (.not.valid) return
        end if
        if (present(transport)) then
          faces=[origin(1,node)>0,origin(1,node)+int(d(1),int64)<domain(1,level+1), &
            origin(2,node)>0,origin(2,node)+int(d(2),int64)<domain(2,level+1), &
            origin(3,node)>0,origin(3,node)+int(d(3),int64)<domain(3,level+1)]
          if (present(inflow_mach)) faces=faces.or.(inflow_mach>=0)
          patch_faces='outflow'
          if(present(boundary_faces)) then
            where(.not.faces.and.(boundary_faces=='slip_wall'.or.boundary_faces=='no_slip_wall'.or. &
              boundary_faces=='isothermal_wall')) &
              patch_faces=boundary_faces
            faces=.true.
          end if
          call advance_reactive_eb_transport_patch_euler_3d(species,transport,snapshot%nodes(node)%geometry, &
            interval,viscosity,conduction,diffusion,barodiffusion,candidate%nodes(node)%state, &
            candidate%nodes(node)%temperature,flux(node),stage_theta,valid,ghost,gt,faces, &
            physical_faces=patch_faces,wall_temperature=wall_temperature)
          if (.not.valid) return
          theta=min(theta,stage_theta)
        else if(present(extended_geometry)) then
          call advance_reactive_eb_plm_ctu_3d(species,working,ghost,gt,vf,interval,flat,candidate%nodes(node)%state, &
            candidate%nodes(node)%temperature,valid,difmag=coefficient,integrated_flux=flux(node),boundary=physical_boundary)
          if(.not.valid) return
          do component=1,nv
            where(snapshot%nodes(node)%geometry%volume_fraction<=0) &
              candidate%nodes(node)%state(component,:,:,:)=snapshot%nodes(node)%state(component,:,:,:)
          end do
          where(snapshot%nodes(node)%geometry%volume_fraction<=0) &
            candidate%nodes(node)%temperature=snapshot%nodes(node)%temperature
          deallocate(flat)
        else
        call advance_reactive_eb_euler_3d(species,snapshot%nodes(node)%state,snapshot%nodes(node)%temperature, &
          snapshot%nodes(node)%geometry,solver,interval,candidate%nodes(node)%state,candidate%nodes(node)%temperature, &
          valid,ghost,gt,flux(node),difmag=coefficient,ghost_volume_fraction=vf)
        end if
        if (.not.valid) return
        deallocate(ghost,gt)
        if(allocated(vf)) deallocate(vf)
        if (present(transport)) cycle
        allocate(source(nv,d(1),d(2),d(3)),exposed(d(1),d(2),d(3)));exposed=.true.
        do child=node+1,size(tree%nodes)
          if (tree%topology%nodes(child)%parent/=node) cycle
          call bounds(tree%topology,child,lo,hi)
          exposed(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
        end do
        call reactive_eb_slip_wall_source_3d(species,snapshot%nodes(node)%state,snapshot%nodes(node)%temperature, &
          snapshot%nodes(node)%geometry,source,valid)
        if (.not.valid) return
        associate(g=>snapshot%nodes(node)%geometry)
        do component=1,nv
          wall(component)=wall(component)+interval*g%dx*g%dy*g%dz* &
            sum(source(component,:,:,:)*g%volume_fraction,mask=exposed)
        end do
        end associate
        deallocate(source,exposed)
      end do
      do node=1,size(tree%nodes)
      do peer=node+1,size(tree%nodes)
        if (tree%topology%nodes(node)%level/=level.or.tree%topology%nodes(peer)%level/=level) cycle
        owner=node;receiver=peer
        do axis=1,3
          if (origin(axis,node)==origin(axis,peer)) cycle
          if (origin(axis,node)>origin(axis,peer)) then
            owner=peer;receiver=node
          end if
          exit
        end do
        offset=origin(:,receiver)-origin(:,owner)
        if (any(offset < -int(tree%topology%nodes(receiver)%shape,int64)).or. &
          any(offset>int(tree%topology%nodes(owner)%shape,int64))) cycle
        call reconcile_reactive_eb_peer_flux_3d(species,tree%nodes(owner)%geometry,tree%nodes(receiver)%geometry, &
          int(offset),flux(owner),flux(receiver),candidate%nodes(receiver)%state, &
          candidate%nodes(receiver)%temperature,valid)
        if (.not.valid) return
      end do
      end do
      ! Reconcile periodic provisional endpoints before descendants interpolate
      ! in time; final fine/coarse synchronization follows the recursive steps.
      if(present(boundary_faces)) then
        associate(g=>tree%nodes(1)%geometry)
        length=[g%x_upper-g%x_lower,g%y_upper-g%y_lower,g%z_upper-g%z_lower]
        end associate
        do seam=1,3
          if(boundary_faces(2*seam-1)/='periodic') cycle
          period=0;period(seam)=length(seam)
          do owner=1,size(tree%nodes)
            if(tree%topology%nodes(owner)%level/=level) cycle
            if(origin(seam,owner)+tree%topology%nodes(owner)%shape(seam)/=domain(seam,level+1)) cycle
            owned=flux(owner)
            do receiver=1,size(tree%nodes)
              if(tree%topology%nodes(receiver)%level/=level.or.origin(seam,receiver)/=0) cycle
              offset=origin(:,receiver)-origin(:,owner)
              offset(seam)=offset(seam)+domain(seam,level+1)
              overlap=.true.
              do tangent=1,3
                if(tangent==seam) cycle
                if(offset(tangent)>=tree%topology%nodes(owner)%shape(tangent).or. &
                  offset(tangent)+tree%topology%nodes(receiver)%shape(tangent)<=0) overlap=.false.
              end do
              if(.not.overlap) cycle
              call reconcile_reactive_eb_peer_flux_3d(species,tree%nodes(owner)%geometry, &
                tree%nodes(receiver)%geometry,int(offset),owned,flux(receiver),candidate%nodes(receiver)%state, &
                candidate%nodes(receiver)%temperature,valid,periodic_length=period)
              if(.not.valid) return
            end do
          end do
        end do
      end if
      do node=1,size(tree%nodes)
        if (tree%topology%nodes(node)%level==level) ends%nodes(node)=candidate%nodes(node)
      end do
      if (level+1<levels) then
        elapsed=0._dp
        substeps=ratios(level+2)
        if (present(transport)) substeps=substeps**2
        do substep=1,substeps
          sub_dt=interval/real(substeps,dp)
          if (substep==substeps) sub_dt=interval-elapsed
          if (sub_dt<=0._dp.or.time+elapsed+sub_dt<=time+elapsed) return
          call advance_level(level+1,time+elapsed,sub_dt,fine_flux,valid)
          if (.not.valid) return
          do node=1,size(tree%nodes)
            if (tree%topology%nodes(node)%level/=level+1) cycle
            flux(node)%x=flux(node)%x+fine_flux(node)%x
            flux(node)%y=flux(node)%y+fine_flux(node)%y
            flux(node)%z=flux(node)%z+fine_flux(node)%z
          end do
          elapsed=elapsed+sub_dt
        end do
      end if
      call synchronize_reactive_eb_tree_flux_3d(species,candidate,flux,valid,level,min(level+1,levels-1),boundary_faces)
      if (.not.valid) return
      good=.true.
    end subroutine advance_level
  end subroutine advance_reactive_eb_tree_subcycle_operator_3d

  ! Serial counterpart of the owner-local adjacent-level repair. Donor
  ! inventory and all internal path faces commit with the fine state; strict
  ! synchronization subsequently applies the matching coarse replacement.
  subroutine repair_reactive_eb_tree_interface_budget_3d(species,solution,flux,ok, &
      minimum_level,boundary_faces,frame,rereflux_delta,rereflux_first_level)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_reactive_tree_3d), intent(inout) :: solution
    type(eb_integrated_flux_3d), intent(inout) :: flux(:)
    integer, intent(in) :: minimum_level
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d), optional, intent(in) :: rereflux_delta
    integer, optional, intent(in) :: rereflux_first_level
    integer :: register_first_level
    ok=.false.
    register_first_level=minimum_level
    if(present(rereflux_first_level)) register_first_level=rereflux_first_level
    if(register_first_level<0.or.register_first_level>minimum_level) return
    call repair_reactive_eb_tree_interface_components_3d(species,solution,flux,ok,minimum_level, &
      boundary_faces,frame,rereflux_delta,0,register_first_level)
  end subroutine repair_reactive_eb_tree_interface_budget_3d

  ! The held scalar moves with every complete gas donor packet, without being
  ! treated as a nonnegative species or reconstructed from total energy.
  subroutine repair_reactive_eb_tree_held_interface_budget_3d(species,solution,internal_energy, &
      flux,internal_energy_flux,ok,minimum_level,boundary_faces,frame,rereflux_delta,internal_energy_rereflux_delta, &
      rereflux_first_level,reactor_energy_bias,native_projection)
    procedure(reactive_eb_funding_admission), optional :: native_projection
    type(amr_tree_state_3d), optional, intent(in) :: reactor_energy_bias
    integer, optional, intent(in) :: rereflux_first_level
    integer :: register_first_level
    type(nasa7_species), intent(in) :: species(:)
    type(eb_reactive_tree_3d), intent(inout) :: solution
    type(amr_tree_state_3d), intent(inout) :: internal_energy
    type(eb_integrated_flux_3d), intent(inout) :: flux(:),internal_energy_flux(:)
    integer, intent(in) :: minimum_level
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d), optional, intent(in) :: rereflux_delta,internal_energy_rereflux_delta
    type(eb_reactive_tree_3d) :: packed,gas_work
    type(amr_tree_state_3d) :: held_work
    type(amr_tree_state_3d), allocatable :: packed_delta
    type(eb_integrated_flux_3d), allocatable :: packets(:),gas_faces(:),held_faces(:)
    integer :: nv,p,n(3),ios
    logical :: good
    ok=.false.;nv=reactive_nvar(size(species))
    if(present(rereflux_delta).neqv.present(internal_energy_rereflux_delta)) return
    register_first_level=minimum_level
    if(present(rereflux_first_level)) register_first_level=rereflux_first_level
    if(register_first_level<0.or.register_first_level>minimum_level) return
    if(.not.valid_tree(solution,nv)) return
    if(.not.internal_energy%is_valid(solution%topology)) return
    if(present(reactor_energy_bias)) then
      if(.not.reactor_energy_bias%is_valid(solution%topology)) return
      do p=1,size(solution%nodes)
        if(size(reactor_energy_bias%fields(p)%values,1)/=1) return
        if(any(.not.ieee_is_finite(reactor_energy_bias%fields(p)%values))) return
      enddo
    endif
    if(size(flux)/=size(solution%nodes).or.size(internal_energy_flux)/=size(solution%nodes)) return
    if(present(rereflux_delta)) then
      if(.not.rereflux_delta%is_valid(solution%topology)) return
      if(.not.internal_energy_rereflux_delta%is_valid(solution%topology)) return
      allocate(packed_delta,stat=ios)
      if(ios/=0) return
      allocate(packed_delta%fields(size(solution%nodes)),stat=ios)
      if(ios/=0) return
    endif
    packed=solution
    allocate(packets(size(solution%nodes)),stat=ios)
    if(ios/=0) return
    do p=1,size(solution%nodes)
      if(size(internal_energy%fields(p)%values,1)/=1) return
      if(.not.all(ieee_is_finite(solution%nodes(p)%state))) return
      if(.not.all(ieee_is_finite(solution%nodes(p)%temperature))) return
      if(.not.all(ieee_is_finite(internal_energy%fields(p)%values))) return
      if(.not.valid_eb_integrated_flux_3d(flux(p),solution%nodes(p)%geometry,nv)) return
      if(.not.valid_eb_integrated_flux_3d(internal_energy_flux(p),solution%nodes(p)%geometry,1)) return
      n=solution%topology%nodes(p)%shape
      deallocate(packed%nodes(p)%state)
      allocate(packed%nodes(p)%state(nv+1,n(1),n(2),n(3)),stat=ios)
      if(ios/=0) return
      packed%nodes(p)%state(:nv,:,:,:)=solution%nodes(p)%state
      packed%nodes(p)%state(nv+1:nv+1,:,:,:)=internal_energy%fields(p)%values
      if(present(rereflux_delta)) then
        if(size(rereflux_delta%fields(p)%values,1)/=nv) return
        if(size(internal_energy_rereflux_delta%fields(p)%values,1)/=1) return
        if(.not.all(ieee_is_finite(rereflux_delta%fields(p)%values))) return
        if(.not.all(ieee_is_finite(internal_energy_rereflux_delta%fields(p)%values))) return
        allocate(packed_delta%fields(p)%values(nv+1,n(1),n(2),n(3)),stat=ios)
        if(ios/=0) return
        packed_delta%fields(p)%values(:nv,:,:,:)=rereflux_delta%fields(p)%values
        packed_delta%fields(p)%values(nv+1:nv+1,:,:,:)=internal_energy_rereflux_delta%fields(p)%values
      endif
      allocate(packets(p)%x(nv+1,0:n(1),n(2),n(3)),packets(p)%y(nv+1,n(1),0:n(2),n(3)), &
        packets(p)%z(nv+1,n(1),n(2),0:n(3)),stat=ios)
      if(ios/=0) return
      packets(p)%x(:nv,:,:,:)=flux(p)%x;packets(p)%x(nv+1:nv+1,:,:,:)=internal_energy_flux(p)%x
      packets(p)%y(:nv,:,:,:)=flux(p)%y;packets(p)%y(nv+1:nv+1,:,:,:)=internal_energy_flux(p)%y
      packets(p)%z(:nv,:,:,:)=flux(p)%z;packets(p)%z(nv+1:nv+1,:,:,:)=internal_energy_flux(p)%z
    enddo
    call repair_reactive_eb_tree_interface_components_3d(species,packed,packets,good,minimum_level, &
      boundary_faces,frame,packed_delta,1,register_first_level,reactor_energy_bias,native_projection)
    if(.not.good) return
    gas_work=solution;held_work=internal_energy;gas_faces=flux;held_faces=internal_energy_flux
    do p=1,size(solution%nodes)
      gas_work%nodes(p)%state=packed%nodes(p)%state(:nv,:,:,:)
      gas_work%nodes(p)%temperature=packed%nodes(p)%temperature
      held_work%fields(p)%values=packed%nodes(p)%state(nv+1:nv+1,:,:,:)
      gas_faces(p)%x=packets(p)%x(:nv,:,:,:);held_faces(p)%x=packets(p)%x(nv+1:nv+1,:,:,:)
      gas_faces(p)%y=packets(p)%y(:nv,:,:,:);held_faces(p)%y=packets(p)%y(nv+1:nv+1,:,:,:)
      gas_faces(p)%z=packets(p)%z(:nv,:,:,:);held_faces(p)%z=packets(p)%z(nv+1:nv+1,:,:,:)
    enddo
    solution=gas_work;call move_alloc(held_work%fields,internal_energy%fields)
    flux=gas_faces;internal_energy_flux=held_faces;ok=.true.
  end subroutine repair_reactive_eb_tree_held_interface_budget_3d

  subroutine repair_reactive_eb_tree_interface_components_3d(species,solution,flux,ok,minimum_level, &
      boundary_faces,frame,rereflux_delta, &
      passive_components,register_first_level,reactor_energy_bias,native_projection)
    procedure(reactive_eb_funding_admission), optional :: native_projection
    type(amr_tree_state_3d), optional, intent(in) :: reactor_energy_bias
    integer, intent(in) :: passive_components,register_first_level
    integer :: active_components
    type(nasa7_species), intent(in) :: species(:)
    type(eb_reactive_tree_3d), intent(inout) :: solution
    type(eb_integrated_flux_3d), intent(inout) :: flux(:)
    integer, intent(in) :: minimum_level
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_patch_tree_3d) :: tree
    type(eb_reactive_tree_3d) :: private_solution
    type(eb_geometry_3d), allocatable :: geometry(:)
    type(amr_tree_state_3d), optional, intent(in) :: rereflux_delta
    type(amr_tree_state_3d) :: state,temperature,candidate,temps,increment,routed_register
    type(eb_integrated_flux_3d), allocatable :: effective(:),working(:)
    real(dp), allocatable :: fine_state(:,:),volumes(:),transfers(:,:),weights(:),areas(:),old(:),base(:),endpoint(:)
    real(dp), allocatable :: corrected(:,:),dc(:),df(:,:),primitive(:)
    real(dp), allocatable :: donor_state(:,:),donor_volume(:),funded_fine(:,:),funded_donor(:,:),exchange(:,:)
    real(dp), allocatable :: packet_values(:,:),low_packet(:),high_packet(:)
    integer, allocatable :: packet_keys(:,:)
    integer :: packet_count,packet_capacity,addition
    integer, allocatable :: used_donors(:,:),next_used_donors(:,:)
    integer, allocatable :: route_node(:),route_cell(:,:),route_axis(:),route_step(:)
    integer :: route_length
    integer, allocatable :: positions(:,:),ratios(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    logical, allocatable :: changed(:)
    integer :: nv,ios,child,parent,other,axis,side,ratio,lo(3),hi(3),n(3),first(3),last(3)
    integer :: i,j,k,x,y,z,q,nf,capacity,face_index(3),coarse_cell(3),fine_cell(3),fine_face(3),fl(3),fh(3)
    integer :: lower(3),upper(3),p,donor_cell(3),internal_face(3)
    integer :: order(3),search_axis,direction,distance,path,donor_axis,donor_step,r,neighbor(3),cursor(3),mirror(3)
    integer :: cursor_node,neighbor_node,donor_node,selected,ancestor_level,ancestor_start
    real(dp) :: coarse_volume,coarse_area,orientation,theta,recovered,speed,maximum_correction,theta_floor
    real(dp) :: coarse_guess
    real(dp) :: coarse_lower,coarse_limit,required_correction,composite_budget(reactive_nvar(size(species))+passive_components)
    logical :: ready,valid,unsafe,covered,local_ok,found,connected_route,periodic(3)
    logical :: route_improved,route_complete,route_can_reroute,immediate_funding
    ok=.false.
    immediate_funding=present(native_projection)
    if(immediate_funding.and..not.present(rereflux_delta)) return
    if(present(reactor_energy_bias).and.passive_components/=1) return
    if(passive_components<0.or.passive_components>1) return
    active_components=reactive_nvar(size(species));nv=active_components+passive_components
    if(.not.valid_tree(solution,nv)) return
    if(present(rereflux_delta)) then
      call route_eb_tree_rereflux_3d(solution,rereflux_delta,routed_register,nv,register_first_level,minimum_level+1,ready)
      if(.not.ready) return
    end if
    call assemble_reactive_eb_tree_reflux_3d(species,solution,flux,increment,effective,ready, &
      minimum_level,minimum_level+1,boundary_faces,frame, &
        passive_components=passive_components)
    if(.not.ready) return
    if(immediate_funding) then
      ! A native source clock must pay its current gas/UEINT deficit before
      ! reaction. Future coarse/fine reflux is not available material stock.
      do p=1,size(solution%nodes)
        increment%fields(p)%values=0
      enddo
      effective=flux
    endif
    block
      ! A negative total after ordinary reflux is just as unfundable as a
      ! negative total with redistribution-register debt.
      private_solution=solution
      do p=1,size(solution%nodes)
        n=solution%topology%nodes(p)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          call add_eb_tree_delta_cell_3d(private_solution%nodes(p)%state(:,i,j,k), &
            increment%fields(p)%values(:,i,j,k),valid)
          if(valid.and.present(rereflux_delta)) call add_eb_tree_delta_cell_3d( &
            private_solution%nodes(p)%state(:,i,j,k),routed_register%fields(p)%values(:,i,j,k),valid)
          if(.not.valid) return
        end do
        end do
        end do
      end do
      call integrate_reactive_eb_tree_3d(private_solution,composite_budget,valid)
      if(.not.valid) return
      if(composite_budget(1)<=0.or.any(composite_budget(6:active_components)<0)) return
    end block
    tree=solution%topology
    allocate(geometry(size(tree%nodes)),state%fields(size(tree%nodes)),temperature%fields(size(tree%nodes)),stat=ios)
    if(ios/=0) return
    do p=1,size(tree%nodes)
      geometry(p)=solution%nodes(p)%geometry
      state%fields(p)%values=solution%nodes(p)%state
      n=tree%nodes(p)%shape
      allocate(temperature%fields(p)%values(1,n(1),n(2),n(3)),stat=ios)
      if(ios/=0) return
      temperature%fields(p)%values(1,:,:,:)=solution%nodes(p)%temperature
    enddo
    periodic=.false.
    if(present(boundary_faces)) periodic=[boundary_faces(1)=='periodic',boundary_faces(3)=='periodic', &
      boundary_faces(5)=='periodic']
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    if(.not.valid) return
    candidate=state;temps=temperature;working=flux
    allocate(old(nv),base(nv),endpoint(nv),low_packet(nv),high_packet(nv), &
      primitive(reactive_nprim(size(species))),changed(size(tree%nodes)),stat=ios)
    ready=(ios==0)
    if(.not.ready) return
    changed=.false.
    ! Apply the pending peer-face divergence to the private fine baseline before
    ! spending donor inventories. Publishing canonical fine faces prevents the
    ! subsequent strict synchronization from applying that divergence twice.
    do p=1,size(tree%nodes)
      if(tree%nodes(p)%level/=minimum_level+1) cycle
      candidate%fields(p)%values=candidate%fields(p)%values+increment%fields(p)%values
      changed(p)=any(increment%fields(p)%values/=0)
      working(p)=effective(p)
    enddo
    if(present(rereflux_delta)) then
      ! Donor inventories include the pending source register once. Descendant
      ! receivers are routed before applying the same conserved flux packets.
      local_ok=.true.
      do p=1,size(tree%nodes)
        if(tree%nodes(p)%level<=minimum_level) cycle
        n=tree%nodes(p)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          call add_eb_tree_delta_cell_3d(candidate%fields(p)%values(:,i,j,k), &
            routed_register%fields(p)%values(:,i,j,k),valid)
          local_ok=local_ok.and.valid
        end do
        end do
        end do
        changed(p)=changed(p).or.any(routed_register%fields(p)%values/=0)
      end do
      ready=local_ok
      if(.not.ready) return
    end if
    if(.not.immediate_funding) then
    do child=1,size(tree%nodes)
      if(tree%nodes(child)%level/=minimum_level+1) cycle
      parent=tree%nodes(child)%parent
      ratio=tree%nodes(child)%patch%refinement_ratio
      capacity=ratio*ratio;n=tree%nodes(parent)%shape
      call bounds(tree,child,lo,hi)
      allocate(fine_state(nv,capacity),volumes(capacity),transfers(nv,capacity),weights(capacity),areas(capacity), &
        positions(3,capacity), &
        donor_state(nv,capacity),donor_volume(capacity),stat=ios)
      ready=(ios==0)
      if(.not.ready) return
      do axis=1,3
      do side=0,1
        orientation=real(2*side-1,dp)
        first=lo;last=hi
        first(axis)=merge(lo(axis)-1,hi(axis)+1,side==0);last(axis)=first(axis)
        do k=first(3),last(3)
        do j=first(2),last(2)
        do i=first(1),last(1)
          coarse_cell=[i,j,k]
          if(any(coarse_cell<1).or.any(coarse_cell>n)) cycle
          if(geometry(parent)%volume_fraction(i,j,k)<=0) cycle
          covered=.false.
          do other=1,size(tree%nodes)
            if(tree%nodes(other)%parent/=parent) cycle
            call bounds(tree,other,lower,upper)
            if(all(coarse_cell>=lower).and.all(coarse_cell<=upper)) covered=.true.
          enddo
          if(covered) cycle
          face_index=coarse_cell;face_index(axis)=face_index(axis)-side
          coarse_volume=geometry(parent)%volume_fraction(i,j,k)* &
            geometry(parent)%dx*geometry(parent)%dy*geometry(parent)%dz
          coarse_area=open_area(parent,axis,face_index)
          if(coarse_area<=0) cycle
          unsafe=.false.
            call cell_increment(parent,coarse_cell,[0,0,0],0,endpoint)
            endpoint=candidate%fields(parent)%values(:,i,j,k)+endpoint
            unsafe=any(endpoint(6:active_components)<0).or.endpoint(1)<=0
          if(.not.unsafe) cycle
            call cell_increment(parent,coarse_cell,face_index,axis,base)
            base=candidate%fields(parent)%values(:,i,j,k)+base
            old=orientation*coarse_area*face_value(working(parent),axis,face_index)
            if(present(rereflux_delta)) then
              ! Encode the coarse register as an offset of this face's old
              ! extensive transfer. The provisional base remains physical;
              ! base + fine - old still equals the complete signed endpoint.
              old=old-coarse_volume*routed_register%fields(parent)%values(:,i,j,k)
            end if
          fl=(coarse_cell-lo)*ratio+1;fh=fl+ratio-1
          fl(axis)=merge(1,tree%nodes(child)%shape(axis),side==0);fh(axis)=fl(axis)
          nf=0
          do z=fl(3),fh(3)
          do y=fl(2),fh(2)
          do x=fl(1),fh(1)
            fine_cell=[x,y,z];fine_face=fine_cell;fine_face(axis)=fine_face(axis)-(1-side)
            if(geometry(child)%volume_fraction(x,y,z)<=0) cycle
            if(open_area(child,axis,fine_face)<=0) cycle
            nf=nf+1;positions(:,nf)=fine_cell
            volumes(nf)=geometry(child)%volume_fraction(x,y,z)*geometry(child)%dx*geometry(child)%dy*geometry(child)%dz
            areas(nf)=open_area(child,axis,fine_face);weights(nf)=areas(nf)/coarse_area
              fine_state(:,nf)=candidate%fields(child)%values(:,x,y,z)
              transfers(:,nf)=orientation*areas(nf)*face_value(effective(child),axis,fine_face)
          enddo
          enddo
          enddo
          valid=nf>0
          ready=(valid)
          if(.not.ready) return
          theta_floor=0
          do q=1,nf
            maximum_correction=1;low_packet=0;high_packet=0
            if(side==0) low_packet=transfers(:,q)-weights(q)*old
            if(side==1) high_packet=weights(q)*old-transfers(:,q)
            call propagate_cell_packet(child,positions(:,q),axis,low_packet,high_packet,.true.,valid,maximum_correction)
            ready=valid
            if(.not.ready) return
            theta_floor=max(theta_floor,1-maximum_correction)
          enddo
          theta=1
          call reconcile_reactive_eb_interface_budget_3d(base,coarse_volume,old,fine_state(:,:nf),volumes(:nf), &
            transfers(:,:nf),weights(:nf),corrected,dc,df,theta,valid,theta_floor, &
        passive_components=passive_components)
          ready=(valid)
          if(.not.ready) then
            call fund_exposed_coarse_endpoint(valid)
            if(valid) then
              call cell_increment(parent,coarse_cell,face_index,axis,base)
              base=candidate%fields(parent)%values(:,i,j,k)+base
              call reconcile_reactive_eb_interface_budget_3d(base,coarse_volume,old,fine_state(:,:nf),volumes(:nf), &
                transfers(:,:nf),weights(:nf),corrected,dc,df,theta,valid,theta_floor, &
        passive_components=passive_components)
              ready=valid
            endif
          endif
          if(.not.ready) then
            ! Each subface chooses its own connected path. A shared donor
            ! is debited immediately in the private state, so later searches
            ! see its remaining inventory. Any failure rolls back the group
            ! and every earlier group before public arrays are published.
            coarse_lower=0;coarse_limit=1
            call reactive_eb_exchange_interval_3d(reshape(base,[nv,1]), &
              reshape((sum(transfers(:,:nf),dim=2)-old)/coarse_volume,[nv,1]),coarse_lower,coarse_limit,valid, &
        passive_components=passive_components)
            if(.not.valid.or.coarse_lower>0) return
            required_correction=1-coarse_limit*(1-32*epsilon(1._dp))
            order=[axis,mod(axis,3)+1,mod(axis+1,3)+1]
            do selected=1,nf
              if(allocated(used_donors)) deallocate(used_donors)
              allocate(used_donors(4,0),stat=ios)
              ready=ios==0
              if(.not.ready) return
              funding: do
              donor_state(:,:nf)=fine_state(:,:nf);donor_volume(:nf)=volumes(:nf)
              call supply_reactive_eb_interface_donors_3d(base,coarse_volume,old,fine_state(:,:nf),volumes(:nf), &
                transfers(:,:nf),weights(:nf),donor_state(:,:nf),donor_volume(:nf),funded_fine,funded_donor,exchange, &
                valid,target_subface=selected,force_exchange=covering_child(child,positions(:,selected))/=0, &
        passive_components=passive_components)
              ready=valid
              if(covering_child(child,positions(:,selected))/=0) then
                call recipient_allowance(selected,maximum_correction,valid)
                if(.not.valid) return
                ready=ready.and.maximum_correction>=required_correction
              endif
              if(ready) exit funding
              found=.false.;connected_route=.false.
              search: do search_axis=1,3
                donor_axis=order(search_axis)
                do direction=0,1
                  donor_step=2*direction-1
                  if(donor_axis==axis) donor_step=(2*direction-1)*(2*side-1)
                  do distance=1,int(min(domain(donor_axis,minimum_level+2)-1,int(huge(distance)-1,int64)))
                    cursor=positions(:,selected);cursor_node=child;valid=.true.
                    do path=1,distance
                      call path_neighbor(cursor_node,cursor,donor_axis,donor_step,neighbor_node,neighbor, &
                        internal_face,mirror,valid)
                      if(.not.valid) exit
                      cursor=neighbor;cursor_node=neighbor_node
                    enddo
                    if(.not.valid) cycle
                    donor_cell=cursor;donor_node=cursor_node
                    do r=1,nf
                      if(donor_node==child.and.all(donor_cell==positions(:,r))) valid=.false.
                    enddo
                    do r=1,size(used_donors,2)
                      if(all(used_donors(:,r)==[donor_node,donor_cell])) valid=.false.
                    enddo
                    if(.not.valid) cycle
                    x=donor_cell(1);y=donor_cell(2);z=donor_cell(3)
                    donor_volume(selected)=cell_volume(donor_node,donor_cell)
                    donor_state(:,selected)=candidate%fields(donor_node)%values(:,x,y,z)
                    call supply_reactive_eb_interface_donors_3d(base,coarse_volume,old,fine_state(:,:nf),volumes(:nf), &
                      transfers(:,:nf),weights(:nf),donor_state(:,:nf),donor_volume(:nf),funded_fine,funded_donor,exchange, &
                      valid,target_subface=selected,force_exchange=covering_child(child,positions(:,selected))/=0, &
        passive_components=passive_components)
                    ready=valid
                    if(.not.ready) cycle
                    call build_straight_route(selected,distance,donor_axis,donor_step,valid)
                    if(.not.valid) return
                    call route_accepts_packet(selected,valid)
                    if(.not.valid.or..not.route_improved) cycle
                    if(route_complete) call minimize_recipient_packet(selected,valid)
                    if(.not.valid) cycle
                    connected_route=.true.;found=.true.;exit search
                  enddo
                enddo
              enddo search
              if(.not.found) then
                call search_connected_donor(selected,found,valid)
                if(.not.valid) return
                if(.not.found) call search_connected_donor(selected,found,valid,allow_partial=.true.)
                if(.not.valid.or..not.found) return
                connected_route=.true.
              endif
              allocate(next_used_donors(4,size(used_donors,2)+1),stat=ios)
              ready=ios==0
              if(.not.ready) return
              next_used_donors(:,:size(used_donors,2))=used_donors
              next_used_donors(:,size(next_used_donors,2))=[donor_node,donor_cell]
              call move_alloc(next_used_donors,used_donors)
              fine_state(:,:nf)=funded_fine;fine_cell=positions(:,selected)
              x=donor_cell(1);y=donor_cell(2);z=donor_cell(3)
              if(covering_child(donor_node,donor_cell)==0) then
                candidate%fields(donor_node)%values(:,x,y,z)=funded_donor(:,selected)
                changed(donor_node)=.true.
              endif
              if(covering_child(child,fine_cell)==0) then
                candidate%fields(child)%values(:,fine_cell(1),fine_cell(2),fine_cell(3))=funded_fine(:,selected)
                changed(child)=.true.
              endif
              cursor=fine_cell;cursor_node=child
              route_length=distance
              if(connected_route) then
                route_length=size(route_node)
              else
                if(allocated(route_node)) deallocate(route_node,route_cell,route_axis,route_step)
                allocate(route_node(route_length),route_cell(3,route_length),route_axis(route_length), &
                  route_step(route_length),stat=ios)
                ready=ios==0
                if(.not.ready) return
              endif
              do path=1,route_length
                if(connected_route) then
                  cursor_node=route_node(path);cursor=route_cell(:,path)
                  donor_axis=route_axis(path);donor_step=route_step(path)
                else
                  route_node(path)=cursor_node;route_cell(:,path)=cursor
                  route_axis(path)=donor_axis;route_step(path)=donor_step
                endif
                call path_neighbor(cursor_node,cursor,donor_axis,donor_step,neighbor_node,neighbor, &
                  internal_face,mirror,valid)
                if(.not.valid) return
                  call add_path_packet(cursor_node,donor_axis,internal_face,-real(donor_step,dp)*exchange(:,selected))
                if(neighbor_node/=cursor_node.or.any(mirror/=internal_face)) then
                    call add_path_packet(neighbor_node,donor_axis,mirror,-real(donor_step,dp)*exchange(:,selected))
                endif
                cursor=neighbor;cursor_node=neighbor_node
              enddo
              call propagate_route_interior(exchange(:,selected),valid)
              if(.not.valid) return
              enddo funding
            enddo
            theta_floor=0
            do q=1,nf
              maximum_correction=1;low_packet=0;high_packet=0
              if(side==0) low_packet=transfers(:,q)-weights(q)*old
              if(side==1) high_packet=weights(q)*old-transfers(:,q)
              call propagate_cell_packet(child,positions(:,q),axis,low_packet,high_packet,.true.,valid,maximum_correction)
              ready=valid
              if(.not.ready) return
              theta_floor=max(theta_floor,1-maximum_correction)
            enddo
            call reconcile_reactive_eb_interface_budget_3d(base,coarse_volume,old,fine_state(:,:nf),volumes(:nf), &
              transfers(:,:nf),weights(:nf),corrected,dc,df,theta,valid,theta_floor, &
        passive_components=passive_components)
            ready=valid
            if(.not.ready) return
          endif
          if(present(rereflux_delta)) then
            coarse_guess=temperature%fields(parent)%values(1,i,j,k)
            call constrain_coarse_eos(valid)
            if(.not.valid) return
          end if
          if(theta==1) cycle
          packet_capacity=0
          do q=1,nf
            addition=transfer_cell_count(child,positions(:,q))
            if(addition<0.or.addition>(huge(packet_capacity)-packet_capacity)/2) return
            packet_capacity=packet_capacity+2*addition
          enddo
          allocate(packet_keys(5,packet_capacity),packet_values(nv,packet_capacity),stat=ios)
          ready=ios==0
          if(.not.ready) return
          packet_count=0
          do q=1,nf
            fine_cell=positions(:,q);fine_face=fine_cell;fine_face(axis)=fine_face(axis)-(1-side)
            low_packet=0;high_packet=0
            if(side==0) low_packet=volumes(q)*df(:,q)
            if(side==1) high_packet=-volumes(q)*df(:,q)
            call propagate_cell_packet(child,fine_cell,axis,low_packet,high_packet,.true.,valid)
            ready=valid
            if(.not.ready) return
              call set_face(working(child),axis,fine_face,orientation*corrected(:,q)/areas(q))
              call set_face(effective(child),axis,fine_face,face_value(working(child),axis,fine_face))
          enddo
          deallocate(packet_keys,packet_values)
          call set_face(effective(parent),axis,face_index, &
            orientation*sum(corrected,dim=2)/coarse_area)
        enddo
        enddo
        enddo
      enddo
      enddo
      deallocate(fine_state,volumes,transfers,weights,areas,positions,donor_state,donor_volume)
    enddo
    endif
    ! StateRedist debt also reaches exposed cells away from fine interfaces
    ! on the closing level and older ancestors. The registered path scans all
    ! such endpoints; flux-only reconciliation keeps its established range.
    ! Candidate state remains pre-register and strict closure charges the
    ! original read-only debt once after these conservative gas/held packets.
    ancestor_start=minimum_level
    do ancestor_level=ancestor_start,register_first_level,-1
      do parent=1,size(tree%nodes)
        if(tree%nodes(parent)%level/=ancestor_level) cycle
        n=tree%nodes(parent)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          coarse_cell=[i,j,k]
          if(geometry(parent)%volume_fraction(i,j,k)<=0) cycle
          if(covering_child(parent,coarse_cell)/=0) cycle
          call cell_increment(parent,coarse_cell,[0,0,0],0,endpoint)
          endpoint=candidate%fields(parent)%values(:,i,j,k)+endpoint
          unsafe=endpoint(1)<=0.or.any(endpoint(6:active_components)<0)
          ! The new flux-only closing-level scan funds thermal EOS deficits;
          ! existing density/species reconciliation keeps its original range.
          if(.not.present(rereflux_delta).and.ancestor_level==minimum_level.and.unsafe) then
            unsafe=.false.
          else if(.not.unsafe) then
            call recover_repair_cell(species,geometry(parent),coarse_cell,endpoint, &
              temps%fields(parent)%values(1,i,j,k),primitive,recovered,speed,valid,frame, &
          repair_node=parent)
            unsafe=.not.valid
          endif
          if(.not.unsafe) cycle
          coarse_volume=cell_volume(parent,coarse_cell)
          call fund_exposed_coarse_endpoint(valid)
          if(.not.valid) return
        enddo
        enddo
        enddo
      enddo
    enddo
    ! Fine candidate inventories already contain routed register debt. Fund
    ! exposed physical leaves before EOS publication, retaining both pre/post
    ! states so strict closure can still charge the original debt once.
      do parent=1,size(tree%nodes)
        if(tree%nodes(parent)%level<=minimum_level) cycle
        n=tree%nodes(parent)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          coarse_cell=[i,j,k]
          if(cell_volume(parent,coarse_cell)<=0.or.covering_child(parent,coarse_cell)/=0) cycle
          unsafe=.false.
          ! Registered physical source debt can exhaust density or species.
          ! Raw zero-debt thermal funding still requires a physical endpoint.
          if(present(rereflux_delta)) then
            if(any(routed_register%fields(parent)%values(:,i,j,k)/=0)) &
              unsafe=candidate%fields(parent)%values(1,i,j,k)<=0.or. &
                any(candidate%fields(parent)%values(6:active_components,i,j,k)<0)
          endif
          if(candidate%fields(parent)%values(1,i,j,k)>0.and. &
            all(candidate%fields(parent)%values(6:active_components,i,j,k)>=0)) then
            call recover_repair_cell(species,geometry(parent),coarse_cell,candidate%fields(parent)%values(:,i,j,k), &
              temps%fields(parent)%values(1,i,j,k),primitive,recovered,speed,valid,frame, &
          repair_node=parent)
            unsafe=.not.valid
          endif
          if(.not.unsafe) cycle
          coarse_volume=cell_volume(parent,coarse_cell)
          call fund_exposed_fine_endpoint(valid)
          ready=valid
          if(.not.ready) return
        enddo
        enddo
        enddo
      enddo
    ! Verify every exposed coarse endpoint after all coupled face adjustments,
    ! and recover changed fine gas states before publishing anything.
    local_ok=.true.
    do p=1,size(tree%nodes)
      n=tree%nodes(p)%shape
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(geometry(p)%volume_fraction(i,j,k)<=0) cycle
        if(tree%nodes(p)%level>=register_first_level.and.tree%nodes(p)%level<=minimum_level) then
          covered=.false.
          do child=1,size(tree%nodes)
            if(tree%nodes(child)%parent/=p) cycle
            call bounds(tree,child,lo,hi)
            if(all([i,j,k]>=lo).and.all([i,j,k]<=hi)) covered=.true.
          enddo
          if(covered) cycle
          call cell_increment(p,[i,j,k],[0,0,0],0,endpoint)
          endpoint=candidate%fields(p)%values(:,i,j,k)+endpoint
          local_ok=local_ok.and.all(ieee_is_finite(endpoint)).and.endpoint(1)>0.and.all(endpoint(6:active_components)>=0)
          call recover_repair_cell(species,geometry(p),[i,j,k],endpoint, &
            temps%fields(p)%values(1,i,j,k),primitive,recovered,speed,valid,frame, &
          repair_node=p)
          local_ok=local_ok.and.valid
        endif
        if(.not.changed(p)) cycle
        if(any(candidate%fields(p)%values(6:active_components,i,j,k)<0)) then
          local_ok=.false.;cycle
        endif
        call recover_repair_cell(species,geometry(p),[i,j,k],candidate%fields(p)%values(:,i,j,k), &
          temperature%fields(p)%values(1,i,j,k),primitive,recovered,speed,valid,frame, &
          repair_node=p)
        if(.not.valid) then
          local_ok=.false.;cycle
        endif
        temps%fields(p)%values(1,i,j,k)=recovered
      enddo
      enddo
      enddo
    enddo
    ready=(local_ok)
    if(.not.ready) return
    if(present(rereflux_delta)) then
      ! Publish the repaired pre-register fine state. Strict closure will apply
      ! the original register exactly once; it remains unconsumed on failure.
      local_ok=.true.
      do p=1,size(tree%nodes)
        if(tree%nodes(p)%level<=minimum_level) cycle
        n=tree%nodes(p)%shape
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(all(routed_register%fields(p)%values(:,i,j,k)==0)) cycle
          call add_eb_tree_delta_cell_3d(candidate%fields(p)%values(:,i,j,k), &
            -routed_register%fields(p)%values(:,i,j,k),valid)
          if(.not.valid) then
            local_ok=.false.;cycle
          end if
          if(geometry(p)%volume_fraction(i,j,k)<=0) cycle
          if(any(candidate%fields(p)%values(6:active_components,i,j,k)<0)) then
            local_ok=.false.;cycle
          end if
          call recover_repair_cell(species,geometry(p),[i,j,k],candidate%fields(p)%values(:,i,j,k), &
            temps%fields(p)%values(1,i,j,k),primitive,recovered,speed,valid,frame, &
          repair_node=p)
          local_ok=local_ok.and.valid
          if(valid) temps%fields(p)%values(1,i,j,k)=recovered
        end do
        end do
        end do
      end do
      ready=local_ok
      if(.not.ready) return
    end if
    private_solution=solution
    do p=1,size(tree%nodes)
      private_solution%nodes(p)%state=candidate%fields(p)%values
      private_solution%nodes(p)%temperature=temps%fields(p)%values(1,:,:,:)
    enddo
    solution=private_solution;flux=working;ok=.true.
  contains
    subroutine recover_repair_cell(table,metrics,index,values,guess,primitive,recovered,speed,valid,rotation,repair_node)
      type(nasa7_species), intent(in) :: table(:)
      type(eb_geometry_3d), intent(in) :: metrics
      integer, intent(in) :: index(3),repair_node
      real(dp), intent(in) :: values(:),guess
      real(dp), intent(out) :: primitive(:),recovered,speed
      logical, intent(out) :: valid
      type(rotating_frame_3d), optional, intent(in) :: rotation
      real(dp) :: gas(active_components),delta(active_components),qp(reactive_nprim(size(table))),tg,sound
      if(passive_components==0) then
        call recover_reactive_eb_cell_3d(table,metrics,index,values,guess,primitive,recovered,speed,valid,rotation)
      else
        call recover_reactive_eb_cell_3d(table,metrics,index,values(:active_components),guess, &
          primitive,recovered,speed,valid,rotation,internal_energy_density=values(active_components+1))
      endif
      if(.not.valid.or..not.present(reactor_energy_bias)) return
      ! Independent UEINT remains the returned temperature authority. Gas and
      ! reactor energy admission are additional queries with private outputs.
      gas=values(:active_components)
      call recover_reactive_eb_cell_3d(table,metrics,index,gas,guess,qp,tg,sound,valid,rotation)
      if(.not.valid) return
      delta=0
      delta(5)=reactor_energy_bias%fields(repair_node)%values(1,index(1),index(2),index(3))
      call add_eb_tree_delta_cell_3d(gas,delta,valid)
      if(.not.valid) return
      call recover_reactive_eb_cell_3d(table,metrics,index,gas,guess,qp,tg,sound,valid,rotation)
    end subroutine recover_repair_cell
    subroutine supply_endpoint_packet(target,delta,target_volume,donor,donor_volume, &
        funded_target,funded_donor,packet,complete,valid,passive_components,allow_inventory_packet, &
        donor_node,donor_cell,donor_delta,target_guess,donor_guess,prefer_inventory)
      real(dp), intent(in) :: target(nv),delta(nv),target_volume,donor(nv),donor_volume,donor_delta(nv)
      real(dp), intent(in) :: target_guess,donor_guess
      real(dp), intent(inout) :: funded_target(nv),funded_donor(nv),packet(nv)
      integer, intent(in) :: passive_components,donor_node,donor_cell(3)
      logical, intent(in) :: allow_inventory_packet
      logical, optional, intent(in) :: prefer_inventory
      logical, intent(out) :: complete,valid
      real(dp) :: value(nv),qp(reactive_nprim(size(species))),tg,sound
      logical :: eos_ok,inventory_required
      if(present(reactor_energy_bias)) then
        call supply_reactive_eb_gas_eos_donor_3d(species,geometry(parent),coarse_cell,target_guess, &
          geometry(donor_node),donor_cell,donor_guess,target,delta,target_volume,donor,donor_delta,donor_volume, &
          funded_target,funded_donor,packet,complete,valid,frame,allow_raw_target=all(delta==0), &
          reactor_energy_bias=[reactor_energy_bias%fields(parent)%values(1,coarse_cell(1),coarse_cell(2),coarse_cell(3)), &
            reactor_energy_bias%fields(donor_node)%values(1,donor_cell(1),donor_cell(2),donor_cell(3))], &
          native_projection=native_projection,inventory_nodes=[parent,donor_node])
        return
      endif
      inventory_required=present(native_projection)
      if(present(prefer_inventory)) inventory_required=inventory_required.or.prefer_inventory
      if(inventory_required.and.allow_inventory_packet) then
          ! Similar donor/receiver density makes a contrast packet too small
          ! for a material sink. Search the complete inventory interval with
          ! physical gas and actual held-energy queries before taking a loan.
          ! Zero offset adds no reactor-energy constraint in the held case.
          if(passive_components==0) then
            call supply_reactive_eb_gas_eos_donor_3d(species,geometry(parent),coarse_cell,target_guess, &
              geometry(donor_node),donor_cell,donor_guess,target,delta,target_volume,donor,donor_delta,donor_volume, &
              funded_target,funded_donor,packet,complete,valid,frame,allow_raw_target=all(delta==0), &
              native_projection=native_projection,inventory_nodes=[parent,donor_node])
          else
            call supply_reactive_eb_gas_eos_donor_3d(species,geometry(parent),coarse_cell,target_guess, &
              geometry(donor_node),donor_cell,donor_guess,target,delta,target_volume,donor,donor_delta,donor_volume, &
              funded_target,funded_donor,packet,complete,valid,frame,allow_raw_target=all(delta==0), &
              reactor_energy_bias=[0._dp,0._dp], &
              native_projection=native_projection,inventory_nodes=[parent,donor_node])
          endif
          return
      endif
      call supply_reactive_eb_endpoint_donor_3d(target,delta,target_volume,donor,donor_volume, &
        funded_target,funded_donor,packet,complete,valid,passive_components,allow_inventory_packet)
      if(valid.and..not.complete) return
      if(valid) then
        value=funded_target
        call add_eb_tree_delta_cell_3d(value,delta,eos_ok)
        if(eos_ok) call recover_repair_cell(species,geometry(parent),coarse_cell,value, &
          target_guess,qp,tg,sound,eos_ok,frame, &
          repair_node=parent)
        if(eos_ok) return
      endif
      if(.not.allow_inventory_packet) then
        ! Flux-only gas/species deficits retain the established interface
        ! contrast/covered donor route. New caloric funding requires a
        ! physical density/composition endpoint before thermal admission.
        value=target
        call add_eb_tree_delta_cell_3d(value,delta,eos_ok)
        valid=.false.;complete=.false.
        if(.not.eos_ok) return
        if(value(1)<=0.or.any(value(6:active_components)<0)) return
      endif
      if(passive_components==0) then
        call supply_reactive_eb_gas_eos_donor_3d(species,geometry(parent),coarse_cell,target_guess, &
          geometry(donor_node),donor_cell,donor_guess,target,delta,target_volume,donor,donor_delta,donor_volume, &
          funded_target,funded_donor,packet,complete,valid,frame,allow_raw_target=all(delta==0), &
              native_projection=native_projection,inventory_nodes=[parent,donor_node])
      else
        call supply_reactive_eb_held_eos_donor_3d(species,geometry(parent),coarse_cell,target_guess, &
          geometry(donor_node),donor_cell,donor_guess,target,delta,target_volume,donor,donor_delta,donor_volume, &
          funded_target,funded_donor,packet,complete,valid,frame,allow_raw_target=all(delta==0))
      endif
    end subroutine supply_endpoint_packet
    subroutine fund_exposed_fine_endpoint(funded)
      logical, intent(out) :: funded
      real(dp) :: source(nv),delta(nv),donor(nv),donor_delta(nv),target_next(nv),target_post(nv)
      real(dp) :: donor_next(nv),donor_pre(nv),packet(nv),face_packet(nv),mirror_packet(nv)
      real(dp) :: source_temp,donor_temp,recovered_target,recovered_donor,sound,area,mirror_area,volume
      integer :: a,step,distance,h,limit,node,next_node,cell(3),next_cell(3),face(3),mirror(3),terminal,terminal_cell(3)
      logical :: connected,path_ok,packet_ok,complete,good,post_good,faces_ok,covered_terminal
      funded=.false.
      search: do a=1,3
      do step=-1,1,2
        limit=int(min(domain(a,tree%nodes(parent)%level+1)-1,int(huge(distance)-1,int64)))
        do distance=1,limit
          node=parent;cell=coarse_cell;path_ok=.true.
          do h=1,distance
            call path_neighbor(node,cell,a,step,next_node,next_cell,face,mirror,connected)
            if(.not.connected) then
              path_ok=.false.;exit
            endif
            node=next_node;cell=next_cell
          enddo
          if(.not.path_ok) exit
          terminal=node;terminal_cell=cell
          source=candidate%fields(parent)%values(:,i,j,k)
          delta=0
          if(present(rereflux_delta)) delta=routed_register%fields(parent)%values(:,i,j,k)
          source_temp=temps%fields(parent)%values(1,i,j,k)
          call add_eb_tree_delta_cell_3d(source,-delta,good)
          if(.not.good) return
          donor=candidate%fields(terminal)%values(:,cell(1),cell(2),cell(3))
          donor_delta=0
          if(present(rereflux_delta)) donor_delta=routed_register%fields(terminal)%values(:,cell(1),cell(2),cell(3))
          donor_temp=temps%fields(terminal)%values(1,cell(1),cell(2),cell(3))
          covered_terminal=covering_child(terminal,cell)/=0
          if(covered_terminal) then
            donor_pre=donor
            call covered_ancestor_inventory(terminal,cell,donor,good)
            if(.not.good) cycle
            donor_delta=donor
            call add_eb_tree_delta_cell_3d(donor_delta,-donor_pre,good)
            if(.not.good) cycle
          endif
          volume=cell_volume(terminal,cell)
          call supply_endpoint_packet(source,delta,coarse_volume,donor,volume, &
            target_next,donor_next,packet,complete,packet_ok,passive_components,.true., &
            terminal,cell,donor_delta,source_temp,donor_temp)
          if(.not.packet_ok) cycle
          call recover_repair_cell(species,geometry(parent),coarse_cell,target_next, &
            source_temp,primitive,recovered_target,sound,good,frame, &
          repair_node=parent)
          if(.not.good.and.(complete.or.any(delta/=0))) cycle
          target_post=target_next
          call add_eb_tree_delta_cell_3d(target_post,delta,good)
          if(.not.good) cycle
          call recover_repair_cell(species,geometry(parent),coarse_cell,target_post, &
            source_temp,primitive,recovered_target,sound,post_good,frame, &
          repair_node=parent)
          if(complete.and..not.post_good) cycle
          call restore_coarse_donor(terminal,cell,donor_next,donor_delta,donor_temp, &
            donor_pre,recovered_donor,good)
          if(.not.good) cycle
          ! Every route face carries the same extensive packet. Exposed
          ! transit cells receive zero net divergence, including peer seams.
          node=parent;cell=coarse_cell;faces_ok=.true.
          do h=1,distance
            call path_neighbor(node,cell,a,step,next_node,next_cell,face,mirror,connected)
            area=open_area(node,a,face);mirror_area=open_area(next_node,a,mirror)
            if(abs(area-mirror_area)>64*epsilon(1._dp)*max(area,mirror_area)) then
              faces_ok=.false.;exit
            endif
            if(area<1) then
              if(any(abs(packet)>huge(1._dp)*area)) then
                faces_ok=.false.;exit
              endif
            endif
            if(mirror_area<1) then
              if(any(abs(packet)>huge(1._dp)*mirror_area)) then
                faces_ok=.false.;exit
              endif
            endif
            face_packet=-real(step,dp)*packet/area
            mirror_packet=-real(step,dp)*packet/mirror_area
            if(any(.not.ieee_is_finite(face_packet)).or.any(.not.ieee_is_finite(mirror_packet))) then
              faces_ok=.false.;exit
            endif
            if(any(abs(face_value(working(node),a,face))>huge(1._dp)-abs(face_packet)).or. &
              any(abs(face_value(effective(node),a,face))>huge(1._dp)-abs(face_packet))) faces_ok=.false.
            if(node/=next_node.or.any(face/=mirror)) then
              if(any(abs(face_value(working(next_node),a,mirror))>huge(1._dp)-abs(mirror_packet)).or. &
                any(abs(face_value(effective(next_node),a,mirror))>huge(1._dp)-abs(mirror_packet))) faces_ok=.false.
            endif
            node=next_node;cell=next_cell
          enddo
          if(.not.faces_ok) cycle
          ! Allocate covered transit/terminal packets to physical descendants.
          ! This checks every changed leaf's pre/post registered EOS and
          ! rolls back a refused route before trying another physical donor.
          call propagate_coarse_route(packet,a,step,distance,good)
          if(.not.good) cycle
          if(covered_terminal) donor_next=donor_pre
          candidate%fields(parent)%values(:,i,j,k)=target_post
          if(post_good) temps%fields(parent)%values(1,i,j,k)=recovered_target
          changed(parent)=.true.
          candidate%fields(terminal)%values(:,terminal_cell(1),terminal_cell(2),terminal_cell(3))=donor_next
          changed(terminal)=.true.
          node=parent;cell=coarse_cell
          do h=1,distance
            call path_neighbor(node,cell,a,step,next_node,next_cell,face,mirror,connected)
            call add_path_packet(node,a,face,-real(step,dp)*packet)
            if(node/=next_node.or.any(face/=mirror)) then
              call add_path_packet(next_node,a,mirror,-real(step,dp)*packet)
            endif
            node=next_node;cell=next_cell
          enddo
          if(complete) then
            funded=.true.;exit search
          endif
        enddo
      enddo
      enddo search
      if(.not.funded.and.present(rereflux_delta)) call fund_adjacent_coarse_endpoint(funded)
    end subroutine fund_exposed_fine_endpoint
    ! A fine leaf can borrow from the exposed physical coarse cell across
    ! its own CF boundary after same-level stock is exhausted. Both faces
    ! carry one extensive gas/held packet; different face areas do not imply
    ! different inventories. Virtual covered coarse cells remain excluded.
    subroutine fund_adjacent_coarse_endpoint(funded)
      logical, intent(out) :: funded
      real(dp) :: source(nv),delta(nv),donor(nv),donor_delta(nv),target_next(nv),target_post(nv)
      real(dp) :: donor_next(nv),donor_pre(nv),packet(nv),face_packet(nv),mirror_packet(nv)
      real(dp) :: source_temp,donor_temp,target_temp,donor_temp_next,sound,area,mirror_area,volume
      integer :: a,step,peer,donor_node,depth,r,cell(3),face(3),mirror(3)
      integer(int64) :: neighbor_global(3),coarse_global(3)
      logical :: packet_ok,complete,good,post_good,faces_ok

      funded=.false.;depth=tree%nodes(parent)%level
      if(depth<1) return
      r=tree%nodes(parent)%patch%refinement_ratio
      if(r<2) return
      do a=1,3
      do step=-1,1,2
        if(step<0.and.coarse_cell(a)/=1) cycle
        if(step>0.and.coarse_cell(a)/=tree%nodes(parent)%shape(a)) cycle
        neighbor_global=origin(:,parent)+int(coarse_cell-1,int64)
        neighbor_global(a)=neighbor_global(a)+step
        if(periodic(a)) neighbor_global(a)=modulo(neighbor_global(a),domain(a,depth+1))
        if(any(neighbor_global<0).or.any(neighbor_global>=domain(:,depth+1))) cycle
        coarse_global=neighbor_global/int(r,int64);donor_node=0
        do peer=1,size(tree%nodes)
          if(tree%nodes(peer)%level/=depth-1) cycle
          if(any(coarse_global<origin(:,peer)).or. &
            any(coarse_global>=origin(:,peer)+int(tree%nodes(peer)%shape,int64))) cycle
          donor_node=peer;cell=int(coarse_global-origin(:,peer))+1;exit
        enddo
        if(donor_node==0) cycle
        if(covering_child(donor_node,cell)/=0) cycle
        volume=cell_volume(donor_node,cell)
        if(volume<=0) cycle
        face=coarse_cell;face(a)=face(a)+min(step,0)
        mirror=cell;mirror(a)=mirror(a)+min(-step,0)
        area=open_area(parent,a,face);mirror_area=open_area(donor_node,a,mirror)
        if(area<=0.or.mirror_area<=0) cycle
        if(area>mirror_area) then
          if(area-mirror_area>64*epsilon(1._dp)*area) cycle
        endif

        source=candidate%fields(parent)%values(:,coarse_cell(1),coarse_cell(2),coarse_cell(3));delta=0
        if(present(rereflux_delta)) delta=routed_register%fields(parent)%values(:,coarse_cell(1),coarse_cell(2),coarse_cell(3))
        source_temp=temps%fields(parent)%values(1,coarse_cell(1),coarse_cell(2),coarse_cell(3))


        call add_eb_tree_delta_cell_3d(source,-delta,good)
        if(.not.good) return

        donor=candidate%fields(donor_node)%values(:,cell(1),cell(2),cell(3))
        call cell_increment(donor_node,cell,[0,0,0],0,donor_delta)
        donor_temp=temps%fields(donor_node)%values(1,cell(1),cell(2),cell(3))


        ! Donor kernels receive its complete post-debt inventory. They
        ! also retain the pre-debt endpoint by subtracting this readonly debt.
        call add_eb_tree_delta_cell_3d(donor,donor_delta,good)
        if(.not.good) cycle
        if(donor(1)<=0.or.any(donor(6:active_components)<0)) cycle
        call supply_endpoint_packet(source,delta,coarse_volume,donor,volume,target_next,donor_next,packet, &
          complete,packet_ok,passive_components,.true.,donor_node,cell,donor_delta,source_temp,donor_temp, &
          prefer_inventory=.true.)
        if(.not.packet_ok) cycle
        call recover_repair_cell(species,geometry(parent),coarse_cell,target_next,source_temp, &
          primitive,target_temp,sound,good,frame,repair_node=parent)
        if(.not.good.and.(complete.or.any(delta/=0))) cycle
        target_post=target_next
        call add_eb_tree_delta_cell_3d(target_post,delta,good)
        if(.not.good) cycle
        call recover_repair_cell(species,geometry(parent),coarse_cell,target_post,source_temp, &
          primitive,target_temp,sound,post_good,frame,repair_node=parent)
        if(complete.and..not.post_good) cycle
        call restore_coarse_donor(donor_node,cell,donor_next,donor_delta,donor_temp,donor_pre,donor_temp_next,good)
        if(.not.good) cycle
        if(area<1) then
          if(any(abs(packet)>huge(1._dp)*area)) cycle
        endif
        if(mirror_area<1) then
          if(any(abs(packet)>huge(1._dp)*mirror_area)) cycle
        endif
        face_packet=-real(step,dp)*packet/area;mirror_packet=-real(step,dp)*packet/mirror_area
        faces_ok=all(ieee_is_finite(face_packet)).and.all(ieee_is_finite(mirror_packet))

        faces_ok=faces_ok.and.all(abs(face_value(working(parent),a,face))<=huge(1._dp)-abs(face_packet)).and. &
          all(abs(face_value(effective(parent),a,face))<=huge(1._dp)-abs(face_packet))


        faces_ok=faces_ok.and.all(abs(face_value(working(donor_node),a,mirror))<=huge(1._dp)-abs(mirror_packet)).and. &
          all(abs(face_value(effective(donor_node),a,mirror))<=huge(1._dp)-abs(mirror_packet))


        if(.not.faces_ok) cycle

        candidate%fields(parent)%values(:,coarse_cell(1),coarse_cell(2),coarse_cell(3))=target_post
        if(post_good) temps%fields(parent)%values(1,coarse_cell(1),coarse_cell(2),coarse_cell(3))=target_temp
        changed(parent)=.true.
        call add_path_packet(parent,a,face,-real(step,dp)*packet)


        candidate%fields(donor_node)%values(:,cell(1),cell(2),cell(3))=donor_pre
        temps%fields(donor_node)%values(1,cell(1),cell(2),cell(3))=donor_temp_next
        changed(donor_node)=.true.
        call add_path_packet(donor_node,a,mirror,-real(step,dp)*packet)

        if(complete) then
          funded=.true.;return
        endif
      enddo
      enddo
    end subroutine fund_adjacent_coarse_endpoint

    subroutine fund_exposed_coarse_endpoint(funded)
      logical, intent(out) :: funded
      real(dp) :: correction(nv),donor_correction(nv),donor_budget(nv),target_next(nv),donor_next(nv),donor_pre(nv),packet(nv)
      real(dp) :: replacement(nv),replacement_effective(nv),other_volume,area,sound,face_packet(nv)
      real(dp) :: mirror_area,mirror_replacement(nv),mirror_effective(nv),mirror_packet(nv)
      real(dp) :: recovered_target,recovered_donor
      integer :: a,step,neighbor(3),face(3),mirror(3),peer
      logical :: complete,packet_complete,valid_target,valid_donor,valid_packet,connected,distant_funded
      logical :: use_inventory
      funded=.false.;complete=.false.;use_inventory=.false.
      if(present(rereflux_delta)) use_inventory=any(routed_register%fields(parent)%values(:,i,j,k)/=0)
      call cell_increment(parent,coarse_cell,[0,0,0],0,correction)
      axes: do a=1,3
      do step=-1,1,2
        neighbor=coarse_cell;neighbor(a)=neighbor(a)+step
        if(any(neighbor<1).or.any(neighbor>n)) cycle
        if(covering_child(parent,neighbor)/=0) cycle
        if(geometry(parent)%volume_fraction(neighbor(1),neighbor(2),neighbor(3))<=0) cycle
        face=coarse_cell;face(a)=face(a)+min(step,0)
        area=open_area(parent,a,face)
        if(area<=0) cycle
        call cell_increment(parent,neighbor,[0,0,0],0,donor_correction)
        call pending_coarse_donor(candidate%fields(parent)%values(:,neighbor(1),neighbor(2),neighbor(3)), &
          donor_correction,donor_budget,valid_donor)
        if(.not.valid_donor) cycle
        other_volume=cell_volume(parent,neighbor)
        target_next=0;donor_next=0;packet=0
        call supply_endpoint_packet(candidate%fields(parent)%values(:,i,j,k),correction, &
          coarse_volume,donor_budget,other_volume, &
          target_next,donor_next,packet,packet_complete,valid_packet, &
        passive_components=passive_components,allow_inventory_packet=use_inventory, &
          donor_node=parent,donor_cell=neighbor,donor_delta=donor_correction, &
          target_guess=temps%fields(parent)%values(1,i,j,k), &
          donor_guess=temps%fields(parent)%values(1,neighbor(1),neighbor(2),neighbor(3)))
        if(.not.valid_packet) cycle
        call recover_repair_cell(species,geometry(parent),coarse_cell,target_next, &
          temps%fields(parent)%values(1,i,j,k),primitive,recovered_target,sound,valid_target,frame, &
          repair_node=parent)
        if(.not.valid_target.and.(packet_complete.or.any(correction/=0))) cycle
        call restore_coarse_donor(parent,neighbor,donor_next,donor_correction, &
          temps%fields(parent)%values(1,neighbor(1),neighbor(2),neighbor(3)),donor_pre,recovered_donor,valid_donor)
        if(.not.valid_donor) cycle
        if(area<1) then
          if(any(abs(packet)>huge(1._dp)*area)) cycle
        endif
        face_packet=-real(step,dp)*packet/area
        if(any(abs(face_value(working(parent),a,face))>huge(1._dp)-abs(face_packet)).or. &
          any(abs(face_value(effective(parent),a,face))>huge(1._dp)-abs(face_packet))) cycle
        replacement=face_value(working(parent),a,face)+face_packet
        replacement_effective=face_value(effective(parent),a,face)+face_packet
        if(any(.not.ieee_is_finite(face_packet)).or.any(.not.ieee_is_finite(replacement)).or. &
          any(.not.ieee_is_finite(replacement_effective))) cycle
        candidate%fields(parent)%values(:,i,j,k)=target_next
        candidate%fields(parent)%values(:,neighbor(1),neighbor(2),neighbor(3))=donor_pre
        if(valid_target) temps%fields(parent)%values(1,i,j,k)=recovered_target
        temps%fields(parent)%values(1,neighbor(1),neighbor(2),neighbor(3))=recovered_donor
        call set_face(working(parent),a,face,replacement)
        call set_face(effective(parent),a,face,replacement_effective)
        changed(parent)=.true.;funded=.true.;complete=packet_complete
        if(complete) exit axes
      enddo
      enddo axes
      if(complete) return
      ! A level can contain several coarse patches. Debit an exposed peer
      ! through both copies of their shared face, including periodic seams.
      peers: do a=1,3
      do step=-1,1,2
        call path_neighbor(parent,coarse_cell,a,step,peer,neighbor,face,mirror,connected)
        if(.not.connected.or.peer==parent) cycle
        if(covering_child(peer,neighbor)/=0) cycle
        area=open_area(parent,a,face);mirror_area=open_area(peer,a,mirror)
        if(abs(area-mirror_area)>64*epsilon(1._dp)*max(area,mirror_area)) cycle
        call cell_increment(peer,neighbor,[0,0,0],0,donor_correction)
        call pending_coarse_donor(candidate%fields(peer)%values(:,neighbor(1),neighbor(2),neighbor(3)), &
          donor_correction,donor_budget,valid_donor)
        if(.not.valid_donor) cycle
        other_volume=cell_volume(peer,neighbor)
        call cell_increment(parent,coarse_cell,[0,0,0],0,correction)
        call supply_endpoint_packet(candidate%fields(parent)%values(:,i,j,k),correction, &
          coarse_volume,donor_budget,other_volume, &
          target_next,donor_next,packet,packet_complete,valid_packet, &
        passive_components=passive_components,allow_inventory_packet=use_inventory, &
          donor_node=peer,donor_cell=neighbor,donor_delta=donor_correction, &
          target_guess=temps%fields(parent)%values(1,i,j,k), &
          donor_guess=temps%fields(peer)%values(1,neighbor(1),neighbor(2),neighbor(3)))
        if(.not.valid_packet) cycle
        call recover_repair_cell(species,geometry(parent),coarse_cell,target_next, &
          temps%fields(parent)%values(1,i,j,k),primitive,recovered_target,sound,valid_target,frame, &
          repair_node=parent)
        if(.not.valid_target.and.(packet_complete.or.any(correction/=0))) cycle
        call restore_coarse_donor(peer,neighbor,donor_next,donor_correction, &
          temps%fields(peer)%values(1,neighbor(1),neighbor(2),neighbor(3)),donor_pre,recovered_donor,valid_donor)
        if(.not.valid_donor) cycle
        if(area<1) then
          if(any(abs(packet)>huge(1._dp)*area)) cycle
        endif
        face_packet=-real(step,dp)*packet/area
        mirror_packet=-real(step,dp)*packet/mirror_area
        if(any(.not.ieee_is_finite(face_packet)).or.any(.not.ieee_is_finite(mirror_packet))) cycle
        if(any(abs(face_value(working(parent),a,face))>huge(1._dp)-abs(face_packet)).or. &
          any(abs(face_value(effective(parent),a,face))>huge(1._dp)-abs(face_packet)).or. &
          any(abs(face_value(working(peer),a,mirror))>huge(1._dp)-abs(mirror_packet)).or. &
          any(abs(face_value(effective(peer),a,mirror))>huge(1._dp)-abs(mirror_packet))) cycle
        replacement=face_value(working(parent),a,face)+face_packet
        replacement_effective=face_value(effective(parent),a,face)+face_packet
        mirror_replacement=face_value(working(peer),a,mirror)+mirror_packet
        mirror_effective=face_value(effective(peer),a,mirror)+mirror_packet
        if(any(.not.ieee_is_finite(replacement)).or.any(.not.ieee_is_finite(replacement_effective)).or. &
          any(.not.ieee_is_finite(mirror_replacement)).or.any(.not.ieee_is_finite(mirror_effective))) cycle
        candidate%fields(parent)%values(:,i,j,k)=target_next
        candidate%fields(peer)%values(:,neighbor(1),neighbor(2),neighbor(3))=donor_pre
        if(valid_target) temps%fields(parent)%values(1,i,j,k)=recovered_target
        temps%fields(peer)%values(1,neighbor(1),neighbor(2),neighbor(3))=recovered_donor
        call set_face(working(parent),a,face,replacement)
        call set_face(effective(parent),a,face,replacement_effective)
        call set_face(working(peer),a,mirror,mirror_replacement)
        call set_face(effective(peer),a,mirror,mirror_effective)
        changed(parent)=.true.;changed(peer)=.true.;funded=.true.;complete=packet_complete
        if(complete) exit peers
      enddo
      enddo peers
      if(.not.complete) then
        call fund_distant_coarse_endpoint(distant_funded,use_inventory)
        funded=funded.or.distant_funded
      endif
    end subroutine fund_exposed_coarse_endpoint
    subroutine fund_distant_coarse_endpoint(funded,use_inventory)
      logical, intent(out) :: funded
      logical, intent(in) :: use_inventory
      real(dp) :: correction(nv),donor_correction(nv),donor_budget(nv),target_next(nv),donor_next(nv),donor_pre(nv),packet(nv)
      real(dp) :: area,mirror_area,volume,face_packet(nv),mirror_packet(nv),sound,recovered_target,recovered_donor
      integer :: a,step,distance,h,limit,node,next_node,cell(3),next_cell(3),face(3),mirror(3)
      integer :: terminal_node,terminal_cell(3)
      logical :: connected,valid_path,valid_packet,complete,valid_target,valid_donor,transit_valid
      funded=.false.
      search: do a=1,3
      do step=-1,1,2
      limit=int(min(domain(a,tree%nodes(parent)%level+1)-1,int(huge(distance)-1,int64)))
      do distance=1,limit
        node=parent;cell=coarse_cell;valid_path=.true.
        do h=1,distance
          call path_neighbor(node,cell,a,step,next_node,next_cell,face,mirror,connected)
          if(.not.connected) then
            valid_path=.false.;exit
          endif
          node=next_node;cell=next_cell
        enddo
        if(.not.valid_path) exit
        terminal_node=node;terminal_cell=cell
        if(distance==1.and.covering_child(terminal_node,terminal_cell)==0) cycle
        ! Adjacent-interface covered cells retain the existing fine-level
        ! route. A farther covered donor is reached through this coarse path;
        ! the interface search does not visit an exposed recipient away from it.
        if(tree%nodes(parent)%level==minimum_level.and.terminal_node==parent.and. &
          covering_child(terminal_node,terminal_cell)/=0) then
          ! The registered native gas/H transaction admits the full packet.
          ! Flux-only reconciliation retains its existing finest-level bounds.
          if(.not.immediate_funding) cycle
        endif
        call cell_increment(parent,coarse_cell,[0,0,0],0,correction)
        call cell_increment(terminal_node,terminal_cell,[0,0,0],0,donor_correction)
        if((tree%nodes(parent)%level<minimum_level.or.present(native_projection)).and. &
          covering_child(terminal_node,terminal_cell)/=0) then
          ! Spend physical descendant inventory after its pending debt once,
          ! including a distant donor on the closing level. Its virtual coarse
          ! average is retained only as the pre-transfer publication anchor.
          call covered_ancestor_inventory(terminal_node,terminal_cell,donor_budget,valid_donor)
          if(valid_donor) then
            donor_correction=donor_budget
            call add_eb_tree_delta_cell_3d(donor_correction, &
              -candidate%fields(terminal_node)%values(:,cell(1),cell(2),cell(3)),valid_donor)
          endif
        else
          call pending_coarse_donor(candidate%fields(terminal_node)%values(:,cell(1),cell(2),cell(3)), &
            donor_correction,donor_budget,valid_donor)
        endif
        if(.not.valid_donor) cycle
        volume=cell_volume(terminal_node,terminal_cell)
        call supply_endpoint_packet(candidate%fields(parent)%values(:,i,j,k),correction, &
          coarse_volume,donor_budget,volume, &
          target_next,donor_next,packet,complete,valid_packet, &
        passive_components=passive_components,allow_inventory_packet=use_inventory, &
          donor_node=terminal_node,donor_cell=terminal_cell,donor_delta=donor_correction, &
          target_guess=temps%fields(parent)%values(1,i,j,k), &
          donor_guess=temps%fields(terminal_node)%values(1,cell(1),cell(2),cell(3)))
        if(.not.valid_packet) cycle
        call recover_repair_cell(species,geometry(parent),coarse_cell,target_next, &
          temps%fields(parent)%values(1,i,j,k),primitive,recovered_target,sound,valid_target,frame, &
          repair_node=parent)
        if(.not.valid_target.and.(complete.or.any(correction/=0))) cycle
        call restore_coarse_donor(terminal_node,terminal_cell,donor_next,donor_correction, &
          temps%fields(terminal_node)%values(1,cell(1),cell(2),cell(3)),donor_pre,recovered_donor,valid_donor)
        if(.not.valid_donor) cycle
        node=parent;cell=coarse_cell;valid_path=.true.
        do h=1,distance
          call path_neighbor(node,cell,a,step,next_node,next_cell,face,mirror,connected)
          area=open_area(node,a,face);mirror_area=open_area(next_node,a,mirror)
          if(abs(area-mirror_area)>64*epsilon(1._dp)*max(area,mirror_area)) then
            valid_path=.false.;exit
          endif
          if(area<1) then
            if(any(abs(packet)>huge(1._dp)*area)) then
              valid_path=.false.;exit
            endif
          endif
          face_packet=-real(step,dp)*packet/area
          mirror_packet=-real(step,dp)*packet/mirror_area
          if(any(.not.ieee_is_finite(face_packet)).or.any(.not.ieee_is_finite(mirror_packet))) then
            valid_path=.false.;exit
          endif
          if(any(.not.ieee_is_finite(face_value(working(node),a,face))).or. &
            any(.not.ieee_is_finite(face_value(effective(node),a,face)))) then
            valid_path=.false.;exit
          endif
          if(any(abs(face_value(working(node),a,face))>huge(1._dp)-abs(face_packet)).or. &
            any(abs(face_value(effective(node),a,face))>huge(1._dp)-abs(face_packet))) then
            valid_path=.false.;exit
          endif
          if(node/=next_node.or.any(face/=mirror)) then
            if(any(.not.ieee_is_finite(face_value(working(next_node),a,mirror))).or. &
              any(.not.ieee_is_finite(face_value(effective(next_node),a,mirror)))) then
              valid_path=.false.;exit
            endif
            if(any(abs(face_value(working(next_node),a,mirror))>huge(1._dp)-abs(mirror_packet)).or. &
              any(abs(face_value(effective(next_node),a,mirror))>huge(1._dp)-abs(mirror_packet))) then
              valid_path=.false.;exit
            endif
          endif
          node=next_node;cell=next_cell
        enddo
        if(.not.valid_path) cycle
        call propagate_coarse_route(packet,a,step,distance,transit_valid)
        if(.not.transit_valid) cycle
        candidate%fields(parent)%values(:,i,j,k)=target_next
        candidate%fields(terminal_node)%values(:,cell(1),cell(2),cell(3))=donor_pre
        if(valid_target) temps%fields(parent)%values(1,i,j,k)=recovered_target
        temps%fields(terminal_node)%values(1,cell(1),cell(2),cell(3))=recovered_donor
        node=parent;cell=coarse_cell
        do h=1,distance
          call path_neighbor(node,cell,a,step,next_node,next_cell,face,mirror,connected)
          call add_path_packet(node,a,face,-real(step,dp)*packet)
          changed(node)=.true.
          if(node/=next_node.or.any(face/=mirror)) then
            call add_path_packet(next_node,a,mirror,-real(step,dp)*packet)
            changed(next_node)=.true.
          endif
          node=next_node;cell=next_cell
        enddo
        funded=.true.
        if(complete) exit search
      enddo
      enddo
      enddo search
    end subroutine fund_distant_coarse_endpoint
    ! Covered ancestor averages are not spendable inventories. Integrate
    ! physical leaves after their pending register/reflux debt, once; fine
    ! leaves already carry that debt privately in this transaction.
    recursive subroutine covered_ancestor_inventory(node,cell,budget,valid)
      integer, intent(in) :: node,cell(3)
      real(dp), intent(out) :: budget(nv)
      logical, intent(out) :: valid
      real(dp) :: child_budget(nv),term(nv),delta(nv),volume,other_volume,total_volume,weight
      real(dp) :: qp(reactive_nprim(size(species))),tg,sound
      integer :: covered_by,lower(3),upper(3),first(3),last(3),r,ii,jj,kk
      budget=0;valid=.false.;volume=cell_volume(node,cell)
      if(volume<=0.or..not.ieee_is_finite(volume)) return
      covered_by=covering_child(node,cell)
      if(covered_by==0) then
        budget=candidate%fields(node)%values(:,cell(1),cell(2),cell(3))
        if(tree%nodes(node)%level<=minimum_level) then
          call cell_increment(node,cell,[0,0,0],0,delta)
          call add_eb_tree_delta_cell_3d(budget,delta,valid)
          if(.not.valid) return
        endif
        valid=all(ieee_is_finite(budget)).and.budget(1)>0.and.all(budget(6:active_components)>=0)
        if(valid) call recover_repair_cell(species,geometry(node),cell,budget, &
          temps%fields(node)%values(1,cell(1),cell(2),cell(3)),qp,tg,sound,valid,frame, &
          repair_node=node)
        return
      endif
      call bounds(tree,covered_by,lower,upper)
      r=tree%nodes(covered_by)%patch%refinement_ratio
      first=(cell-lower)*r+1;last=first+r-1;total_volume=0
      do kk=first(3),last(3)
      do jj=first(2),last(2)
      do ii=first(1),last(1)
        other_volume=cell_volume(covered_by,[ii,jj,kk])
        if(other_volume<=0) cycle
        if(.not.ieee_is_finite(other_volume).or.total_volume>huge(1._dp)-other_volume) then
          valid=.false.;return
        endif
        call covered_ancestor_inventory(covered_by,[ii,jj,kk],child_budget,valid)
        if(.not.valid) return
        total_volume=total_volume+other_volume;weight=other_volume/volume
        if(weight>1) then
          if(any(abs(child_budget)>huge(1._dp)/weight)) then
            valid=.false.;return
          endif
        endif
        term=weight*child_budget
        call add_eb_tree_delta_cell_3d(budget,term,valid)
        if(.not.valid) return
      enddo
      enddo
      enddo
      valid=abs(total_volume-volume)<=64*epsilon(1._dp)*volume
    end subroutine covered_ancestor_inventory
    subroutine qualify_ancestor_route(baseline,valid)
      type(amr_tree_state_3d), intent(in) :: baseline
      logical, intent(out) :: valid
      real(dp) :: value(nv),delta(nv),qp(reactive_nprim(size(species))),tg,sound
      integer :: node,d(3),ii,jj,kk
      logical :: good,local_ok
      local_ok=.true.
      do node=1,size(tree%nodes)
        if(.not.changed(node)) cycle
        d=tree%nodes(node)%shape
        do kk=1,d(3)
        do jj=1,d(2)
        do ii=1,d(1)
          if(cell_volume(node,[ii,jj,kk])<=0.or.covering_child(node,[ii,jj,kk])/=0) cycle
          value=candidate%fields(node)%values(:,ii,jj,kk)
          if(all(value==baseline%fields(node)%values(:,ii,jj,kk))) cycle
          good=all(ieee_is_finite(value)).and.value(1)>0.and.all(value(6:active_components)>=0)
          if(good) call recover_repair_cell(species,geometry(node),[ii,jj,kk],value, &
            temps%fields(node)%values(1,ii,jj,kk),qp,tg,sound,good,frame, &
          repair_node=node)
          local_ok=local_ok.and.good
          if(tree%nodes(node)%level>minimum_level) then
            if(tree%nodes(parent)%level<=minimum_level) cycle
            if(.not.present(rereflux_delta)) cycle
            call add_eb_tree_delta_cell_3d(value,-routed_register%fields(node)%values(:,ii,jj,kk),good)
          else
            call cell_increment(node,[ii,jj,kk],[0,0,0],0,delta)
            call add_eb_tree_delta_cell_3d(value,delta,good)
          endif
          if(good) good=value(1)>0.and.all(value(6:active_components)>=0)
          if(good) call recover_repair_cell(species,geometry(node),[ii,jj,kk],value, &
            temps%fields(node)%values(1,ii,jj,kk),qp,tg,sound,good,frame, &
          repair_node=node)
          local_ok=local_ok.and.good
        enddo
        enddo
        enddo
      enddo
      valid=local_ok
    end subroutine qualify_ancestor_route
    subroutine pending_coarse_donor(raw,delta,budget,valid)
      real(dp), intent(in) :: raw(nv),delta(nv)
      real(dp), intent(out) :: budget(nv)
      logical, intent(out) :: valid
      budget=raw
      call add_eb_tree_delta_cell_3d(budget,delta,valid)
      if(valid) valid=budget(1)>0.and.all(budget(6:active_components)>=0)
    end subroutine pending_coarse_donor
    subroutine restore_coarse_donor(node,cell,post,delta,guess,pre,recovered,valid)
      integer, intent(in) :: node,cell(3)
      real(dp), intent(in) :: post(nv),delta(nv),guess
      real(dp), intent(out) :: pre(nv),recovered
      logical, intent(out) :: valid
      real(dp) :: post_temp,sound
      call recover_repair_cell(species,geometry(node),cell,post,guess,primitive,post_temp,sound,valid,frame, &
          repair_node=node)
      if(.not.valid) return
      pre=post
      call add_eb_tree_delta_cell_3d(pre,-delta,valid)
      if(.not.valid) return
      valid=pre(1)>0.and.all(pre(6:active_components)>=0)
      if(.not.valid) return
      call recover_repair_cell(species,geometry(node),cell,pre,guess,primitive,recovered,sound,valid,frame, &
          repair_node=node)
    end subroutine restore_coarse_donor
    subroutine propagate_coarse_route(packet,a,step,distance,valid)
      real(dp), intent(in) :: packet(nv)
      integer, intent(in) :: a,step,distance
      logical, intent(out) :: valid
      type(amr_tree_state_3d) :: saved_candidate
      type(eb_integrated_flux_3d), allocatable :: saved_working(:),saved_effective(:)
      logical, allocatable :: saved_changed(:)
      real(dp) :: low(nv),high(nv),maximum
      integer :: h,node,next_node,cell(3),next_cell(3),face(3),mirror(3),capacity,count,ios,p
      logical :: connected,trial_ok
      valid=.false.;capacity=0
      node=parent;cell=coarse_cell
      do h=1,distance
        call path_neighbor(node,cell,a,step,next_node,next_cell,face,mirror,connected)
        if(covering_child(next_node,next_cell)/=0) then
          low=-real(step,dp)*packet;high=low
          if(h==distance) then
            low=0;high=0
            if(step>0) low=-packet
            if(step<0) high=packet
          endif
          count=transfer_cell_count(next_node,next_cell)
          if(count<0.or.count>(huge(capacity)-capacity)/6) return
          capacity=capacity+6*count
          maximum=1
          call propagate_cell_packet(next_node,next_cell,a,low,high,.true.,trial_ok,maximum)
          if(.not.trial_ok.or.maximum<1) return
        endif
        node=next_node;cell=next_cell
      enddo
      if(capacity==0) then
        valid=.true.;return
      endif
      allocate(packet_keys(5,capacity),packet_values(nv,capacity),stat=ios)
      if(ios/=0) then
        if(allocated(packet_keys)) deallocate(packet_keys)
        if(allocated(packet_values)) deallocate(packet_values)
        return
      endif
      packet_count=0
      saved_candidate=candidate;saved_working=working;saved_effective=effective;saved_changed=changed
      node=parent;cell=coarse_cell
      do h=1,distance
        call path_neighbor(node,cell,a,step,next_node,next_cell,face,mirror,connected)
        if(covering_child(next_node,next_cell)/=0) then
          low=-real(step,dp)*packet;high=low
          if(h==distance) then
            low=0;high=0
            if(step>0) low=-packet
            if(step<0) high=packet
          endif
          call propagate_cell_packet(next_node,next_cell,a,low,high,.true.,trial_ok)
          if(.not.trial_ok) exit
        endif
        node=next_node;cell=next_cell
      enddo
      valid=trial_ok
      if(valid) then
        do p=1,size(tree%nodes)
          if(.not.changed(p)) cycle
          if(any(.not.ieee_is_finite(candidate%fields(p)%values)).or. &
            any(.not.ieee_is_finite(working(p)%x)).or.any(.not.ieee_is_finite(working(p)%y)).or. &
            any(.not.ieee_is_finite(working(p)%z)).or.any(.not.ieee_is_finite(effective(p)%x)).or. &
            any(.not.ieee_is_finite(effective(p)%y)).or.any(.not.ieee_is_finite(effective(p)%z))) then
            valid=.false.;exit
          endif
        enddo
      endif
      if(valid.and.tree%nodes(parent)%level/=minimum_level) call qualify_ancestor_route(saved_candidate,valid)
      if(.not.valid) then
        candidate=saved_candidate;working=saved_working;effective=saved_effective;changed=saved_changed
      endif
      deallocate(packet_keys,packet_values)
    end subroutine propagate_coarse_route
    ! Tighten the existing conservative face blend when a nearly empty
    ! coarse endpoint cannot satisfy the unchanged mixture EOS contract.
    ! Donor/species bounds stay enforced by the same interval kernel.
    subroutine constrain_coarse_eos(good)
      logical, intent(out) :: good
      real(dp), allocatable :: best_faces(:,:),best_coarse(:),best_fine(:,:)
      real(dp) :: low,high,cap,best_theta
      logical :: accepted,found
      integer :: attempt
      call qualify_coarse_endpoint(good)
      if(good) return
      low=theta_floor;high=theta;found=.false.
      do attempt=1,64
        cap=low+.5_dp*(high-low)
        if(cap<=low.or.cap>=high) exit
        call reconcile_reactive_eb_interface_budget_3d(base,coarse_volume,old,fine_state(:,:nf),volumes(:nf), &
          transfers(:,:nf),weights(:nf),corrected,dc,df,theta,accepted,theta_floor,cap, &
        passive_components=passive_components)
        if(.not.accepted) then
          low=cap
          cycle
        end if
        call qualify_coarse_endpoint(accepted)
        if(accepted) then
          found=.true.;best_theta=theta
          best_faces=corrected;best_coarse=dc;best_fine=df
          low=cap
        else
          high=cap
        end if
      end do
      good=found
      if(.not.good) return
      corrected=best_faces;dc=best_coarse;df=best_fine;theta=best_theta
    end subroutine constrain_coarse_eos
    subroutine qualify_coarse_endpoint(good)
      logical, intent(out) :: good
      type(eb_integrated_flux_3d) :: trial_flux
      real(dp), allocatable :: delta(:,:,:,:)
      logical, allocatable :: refined(:,:,:)
      real(dp) :: trial_value(nv),qp(reactive_nprim(size(species))),tg,sound
      integer :: child,bxlo(3),bxhi(3),d(3),ios
      good=.false.
        d=tree%nodes(parent)%shape
        allocate(refined(d(1),d(2),d(3)),stat=ios)
        if(ios==0) then
          refined=.false.
          do child=1,size(tree%nodes)
            if(tree%nodes(child)%parent/=parent) cycle
            call bounds(tree,child,bxlo,bxhi)
            refined(bxlo(1):bxhi(1),bxlo(2):bxhi(2),bxlo(3):bxhi(3))=.true.
          end do
          trial_flux=effective(parent)
          call set_face(trial_flux,axis,face_index,orientation*sum(corrected,dim=2)/coarse_area)
          call compute_reactive_eb_reflux_increment_3d(geometry(parent),refined,working(parent),trial_flux,delta, &
            good,all_faces=.true.,redistribution_delta=routed_register%fields(parent)%values)
          if(good) then
            trial_value=candidate%fields(parent)%values(:,coarse_cell(1),coarse_cell(2),coarse_cell(3))+ &
              delta(:,coarse_cell(1),coarse_cell(2),coarse_cell(3))
            call recover_repair_cell(species,geometry(parent),coarse_cell,trial_value, &
              coarse_guess,qp,tg,sound,good,frame, &
          repair_node=parent)
          end if
        end if
    end subroutine qualify_coarse_endpoint
    ! Breadth-first fallback over open same-level faces. Geometry and route
    ! metadata are shared; only the donor endpoint vector is broadcast by MPI.
    recursive subroutine search_connected_donor(target,found,valid,allow_partial,excluded_edges,fixed_donor,rejected_route)
      integer, intent(in) :: target
      integer, optional, intent(in) :: excluded_edges(:,:),fixed_donor(4)
      logical, optional, intent(out) :: rejected_route
      logical, intent(out) :: found,valid
      logical, optional, intent(in) :: allow_partial
      integer, allocatable :: queue_node(:),queue_cell(:,:),previous(:),edge_axis(:),edge_step(:)
      integer, allocatable :: new_node(:),new_cell(:,:),new_previous(:),new_axis(:),new_step(:)
      integer(int64), allocatable :: node_offset(:),visited_bits(:)
      integer :: capacity,limit,next_capacity,word,bit,head,tail,entry,node,next_node,cell(3),next_cell(3),face(3),next_face(3)
      integer :: peer,a_index,a,step,direction,r,ios,path_entry
      integer(int64) :: total,product_cells,vertex
      integer :: excluded_peer,excluded_cell(3),excluded_face(3),excluded_mirror(3),blocked
      logical :: edge_ok,protected,ready,alternative_found
      found=.false.;valid=.false.;total=0
      if(present(rejected_route)) rejected_route=.false.
      if(present(excluded_edges)) then
        if(size(excluded_edges,1)/=6) return
      endif
      allocate(node_offset(size(tree%nodes)),stat=ios)
      ready=ios==0
      if(.not.ready) return
      node_offset=-1
      do peer=1,size(tree%nodes)
        if(tree%nodes(peer)%level/=minimum_level+1) cycle
        node_offset(peer)=total;product_cells=1
        do a=1,3
          if(product_cells>int(huge(capacity),int64)/tree%nodes(peer)%shape(a)) return
          product_cells=product_cells*tree%nodes(peer)%shape(a)
        enddo
        if(total>int(huge(capacity),int64)-product_cells) return
        total=total+product_cells
      enddo
      if(total<1) return
      limit=int(total);capacity=min(limit,64)
      ! One bit per possible vertex replaces the scan through every queued
      ! cell. Route arrays grow only when another frontier slot is needed.
      allocate(visited_bits(int((total+63_int64)/64_int64)),queue_node(capacity),queue_cell(3,capacity), &
        previous(capacity),edge_axis(capacity),edge_step(capacity),stat=ios)
      ready=ios==0
      if(.not.ready) return
      visited_bits=0
      head=1;tail=1;queue_node(1)=child;queue_cell(:,1)=positions(:,target);previous(1)=0
      edge_axis(1)=0;edge_step(1)=0;cell=positions(:,target)
      vertex=node_offset(child)+(int(cell(3)-1,int64)*int(tree%nodes(child)%shape(2),int64)+cell(2)-1)* &
        int(tree%nodes(child)%shape(1),int64)+cell(1)-1
      word=int(vertex/64_int64)+1;bit=int(modulo(vertex,64_int64))
      visited_bits(word)=ibset(visited_bits(word),bit)
      do while(head<=tail)
        node=queue_node(head);cell=queue_cell(:,head)
        do a_index=1,3
          a=order(a_index)
          do direction=0,1
            step=2*direction-1
            if(a==axis) step=step*(2*side-1)
            call path_neighbor(node,cell,a,step,next_node,next_cell,face,next_face,edge_ok)
            if(.not.edge_ok) cycle
            if(present(excluded_edges)) then
              protected=.false.
              do blocked=1,size(excluded_edges,2)
                if(a/=excluded_edges(5,blocked)) cycle
                if(node==excluded_edges(1,blocked).and.all(cell==excluded_edges(2:4,blocked)).and. &
                  step==excluded_edges(6,blocked)) then
                  protected=.true.;exit
                endif
                call path_neighbor(excluded_edges(1,blocked),excluded_edges(2:4,blocked), &
                  excluded_edges(5,blocked),excluded_edges(6,blocked),excluded_peer,excluded_cell, &
                  excluded_face,excluded_mirror,edge_ok)
                if(.not.edge_ok) return
                if(node==excluded_peer.and.all(cell==excluded_cell).and.step==-excluded_edges(6,blocked)) then
                  protected=.true.;exit
                endif
              enddo
              if(protected) cycle
            endif
            vertex=node_offset(next_node)+(int(next_cell(3)-1,int64)* &
              int(tree%nodes(next_node)%shape(2),int64)+next_cell(2)-1)* &
              int(tree%nodes(next_node)%shape(1),int64)+next_cell(1)-1
            word=int(vertex/64_int64)+1;bit=int(modulo(vertex,64_int64))
            if(btest(visited_bits(word),bit)) cycle
            if(tail==capacity) then
              if(capacity==limit) return
              next_capacity=limit
              if(capacity<=limit/2) next_capacity=2*capacity
              allocate(new_node(next_capacity),new_cell(3,next_capacity),new_previous(next_capacity), &
                new_axis(next_capacity),new_step(next_capacity),stat=ios)
              ready=ios==0
              if(.not.ready) return
              new_node(:tail)=queue_node;new_cell(:,:tail)=queue_cell;new_previous(:tail)=previous
              new_axis(:tail)=edge_axis;new_step(:tail)=edge_step
              call move_alloc(new_node,queue_node);call move_alloc(new_cell,queue_cell)
              call move_alloc(new_previous,previous);call move_alloc(new_axis,edge_axis)
              call move_alloc(new_step,edge_step);capacity=next_capacity
            endif
            visited_bits(word)=ibset(visited_bits(word),bit)
            tail=tail+1;queue_node(tail)=next_node;queue_cell(:,tail)=next_cell;previous(tail)=head
            edge_axis(tail)=a;edge_step(tail)=step
            protected=.false.
            do r=1,nf
              if(next_node==child.and.all(next_cell==positions(:,r))) protected=.true.
            enddo
            do r=1,size(used_donors,2)
              if(all(used_donors(:,r)==[next_node,next_cell])) protected=.true.
            enddo
            if(protected) cycle
            if(present(fixed_donor)) then
              if(any([next_node,next_cell]/=fixed_donor)) cycle
            endif
            donor_volume(target)=cell_volume(next_node,next_cell)
            donor_state(:,target)=candidate%fields(next_node)%values(:,next_cell(1),next_cell(2),next_cell(3))
            call supply_reactive_eb_interface_donors_3d(base,coarse_volume,old,fine_state(:,:nf),volumes(:nf), &
              transfers(:,:nf),weights(:nf),donor_state(:,:nf),donor_volume(:nf),funded_fine,funded_donor,exchange, &
              edge_ok,target_subface=target,allow_partial=allow_partial, &
              force_exchange=covering_child(child,positions(:,target))/=0, &
        passive_components=passive_components)
            ready=edge_ok
            if(.not.ready) cycle
            donor_node=next_node;donor_cell=next_cell
            if(allocated(route_node)) deallocate(route_node,route_cell,route_axis,route_step)
            route_length=0;path_entry=tail
            do while(previous(path_entry)>0)
              route_length=route_length+1;path_entry=previous(path_entry)
            enddo
            allocate(route_node(route_length),route_cell(3,route_length),route_axis(route_length), &
              route_step(route_length),stat=ios)
            ready=ios==0
            if(.not.ready) return
            path_entry=tail
            do entry=route_length,1,-1
              route_axis(entry)=edge_axis(path_entry);route_step(entry)=edge_step(path_entry)
              path_entry=previous(path_entry)
              route_node(entry)=queue_node(path_entry);route_cell(:,entry)=queue_cell(:,path_entry)
            enddo
            call route_accepts_packet(target,edge_ok,allow_partial)
            if(.not.edge_ok.or..not.route_improved) then
              if(present(rejected_route)) rejected_route=route_can_reroute.and.any(exchange(:,target)/=0)
              if(.not.present(fixed_donor).and.route_can_reroute.and.any(exchange(:,target)/=0)) then
                call find_alternate_donor_route(target,alternative_found,valid,allow_partial)
                if(.not.valid) return
                if(alternative_found) then
                  found=.true.;return
                endif
              endif
              cycle
            endif
            if(route_complete) call minimize_recipient_packet(target,edge_ok)
            if(.not.edge_ok) cycle
            found=.true.;valid=.true.;return
          enddo
        enddo
        head=head+1
      enddo
      valid=.true.
    end subroutine search_connected_donor
    ! Explore canonical exclusion sets in breadth-first order. A rejected
    ! rerouted path contributes its own edges, so a feasible route may avoid
    ! multiple independently inadmissible corridors. Duplicate sets are
    ! visited once; gas/flux publication remains outside this metadata search.
    recursive subroutine find_alternate_donor_route(target,found,valid,allow_partial)
      integer, intent(in) :: target
      logical, intent(out) :: found,valid
      logical, optional, intent(in) :: allow_partial
      type :: retry_state
        integer, allocatable :: blocked(:,:),path(:,:)
      end type retry_state
      type(retry_state), allocatable :: frontier(:),grown(:)
      integer, allocatable :: saved_node(:),saved_cell(:,:),saved_axis(:),saved_step(:),blocked_edges(:,:)
      real(dp) :: saved_fine(nv),saved_donor(nv),saved_exchange(nv),missing(nv),allowance
      integer :: edge,count,ios,endpoint(4),key(6),temporary(6),head,tail,capacity,next_capacity,nblocked,entry,prior
      integer :: incident_node,incident_cell(3),incident_face(3),incident_mirror(3)
      logical :: ready,rejected,duplicate
      found=.false.;valid=.false.;count=route_length;endpoint=[donor_node,donor_cell]
      call recipient_allowance(target,allowance,valid,missing)
      if(.not.valid) return
      ! A covered endpoint receives a nonnegative geometric distribution of
      ! this same packet, independently of the remaining route. If an actual
      ! density/species deficit exists and no deficient component receives a
      ! positive increment, changing incident faces cannot improve its budget.
      ! Positive total mass alone does not replenish an absent species. Keep
      ! the prior conservative gate when no physical deficit is established;
      ! EOS-only rejection is not such a certificate. Held energy is signed
      ! and does not participate in this species-admission certificate.
      if(missing(1)>0.or.any(missing(6:active_components)>0)) then
        if(.not.(missing(1)>0.and.exchange(1,target)>0).and. &
          .not.any(missing(6:active_components)>0.and.exchange(6:active_components,target)>0)) return
      else
        if(exchange(1,target)<=0) return
      endif
      valid=.false.
      if(count>huge(count)/6) return
      capacity=16
      allocate(saved_node(count),saved_cell(3,count),saved_axis(count),saved_step(count),frontier(capacity), &
        stat=ios)
      ready=ios==0
      if(.not.ready) then
        valid=.false.;return
      endif
      saved_node=route_node;saved_cell=route_cell;saved_axis=route_axis;saved_step=route_step
      saved_fine=funded_fine(:,target);saved_donor=funded_donor(:,target);saved_exchange=exchange(:,target)
      allocate(frontier(1)%blocked(6,0),frontier(1)%path(6,count),stat=ios)
      ready=ios==0
      if(.not.ready) then
        valid=.false.;return
      endif
      do edge=1,count
        frontier(1)%path(:,edge)=[saved_node(edge),saved_cell(:,edge),saved_axis(edge),saved_step(edge)]
      enddo
      head=1;tail=1
      do while(head<=tail)
        nblocked=size(frontier(head)%blocked,2)
        if(nblocked>=huge(nblocked)/6) then
          valid=.false.;return
        endif
        do edge=1,size(frontier(head)%path,2)
          ! Rerouting ordinary exposed interiors changes only the face ledger:
          ! their incoming/outgoing complete packets cancel without a state
          ! update. Only edges incident to a covered coarse cell can change
          ! descendant spreading and the admission of this rejected packet.
          ! Retain all such edges, including multiple covered corridors.
          call path_neighbor(frontier(head)%path(1,edge),frontier(head)%path(2:4,edge), &
            frontier(head)%path(5,edge),frontier(head)%path(6,edge),incident_node,incident_cell, &
            incident_face,incident_mirror,valid)
          if(.not.valid) return
          if(covering_child(frontier(head)%path(1,edge),frontier(head)%path(2:4,edge))==0.and. &
            covering_child(incident_node,incident_cell)==0) cycle
          call canonical_route_edge(frontier(head)%path(:,edge),key,valid)
          if(.not.valid) return
          duplicate=.false.
          do entry=1,nblocked
            if(all(key==frontier(head)%blocked(:,entry))) duplicate=.true.
          enddo
          if(duplicate) cycle
          allocate(blocked_edges(6,nblocked+1),stat=ios)
      ready=ios==0
          if(.not.ready) then
        valid=.false.;return
      endif
          blocked_edges(:,:nblocked)=frontier(head)%blocked;blocked_edges(:,nblocked+1)=key
          ! Insertion maintains one representation for every undirected set.
          do entry=nblocked+1,2,-1
            if(.not.route_edge_less(blocked_edges(:,entry),blocked_edges(:,entry-1))) exit
            temporary=blocked_edges(:,entry-1);blocked_edges(:,entry-1)=blocked_edges(:,entry)
            blocked_edges(:,entry)=temporary
          enddo
          duplicate=.false.
          do prior=1,tail
            if(size(frontier(prior)%blocked,2)/=nblocked+1) cycle
            if(all(frontier(prior)%blocked==blocked_edges)) then
              duplicate=.true.;exit
            endif
          enddo
          if(duplicate) then
            deallocate(blocked_edges);cycle
          endif
          if(tail==capacity) then
            if(capacity==huge(capacity)) then
              valid=.false.;return
            endif
            next_capacity=huge(capacity)
            if(capacity<=huge(capacity)/2) next_capacity=2*capacity
            allocate(grown(next_capacity),stat=ios)
      ready=ios==0
            if(.not.ready) then
        valid=.false.;return
      endif
            do prior=1,tail
              call move_alloc(frontier(prior)%blocked,grown(prior)%blocked)
              call move_alloc(frontier(prior)%path,grown(prior)%path)
            enddo
            call move_alloc(grown,frontier);capacity=next_capacity
          endif
          tail=tail+1;call move_alloc(blocked_edges,frontier(tail)%blocked)
          call search_connected_donor(target,found,valid,allow_partial,frontier(tail)%blocked,endpoint,rejected)
          if(.not.valid.or.found) return
          count=0
          if(rejected) count=route_length
          if(count>huge(count)/6) then
            valid=.false.;return
          endif
          allocate(frontier(tail)%path(6,count),stat=ios)
      ready=ios==0
          if(.not.ready) then
        valid=.false.;return
      endif
          do entry=1,count
            frontier(tail)%path(:,entry)=[route_node(entry),route_cell(:,entry),route_axis(entry),route_step(entry)]
          enddo
        enddo
        head=head+1
      enddo
      ! Restore the failed original candidate before resuming other donors.
      count=size(saved_node)
      if(allocated(route_node)) deallocate(route_node,route_cell,route_axis,route_step)
      allocate(route_node(count),route_cell(3,count),route_axis(count),route_step(count),stat=ios)
      ready=ios==0
      if(.not.ready) then
        valid=.false.;return
      endif
      route_node=saved_node;route_cell=saved_cell;route_axis=saved_axis;route_step=saved_step;route_length=count
      donor_node=endpoint(1);donor_cell=endpoint(2:4)
      funded_fine(:,target)=saved_fine;funded_donor(:,target)=saved_donor;exchange(:,target)=saved_exchange
      valid=.true.
    end subroutine find_alternate_donor_route
    subroutine canonical_route_edge(edge,key,valid)
      integer, intent(in) :: edge(6)
      integer, intent(out) :: key(6)
      logical, intent(out) :: valid
      integer :: next_node,next_cell(3),face(3),mirror(3),reverse(6)
      call path_neighbor(edge(1),edge(2:4),edge(5),edge(6),next_node,next_cell,face,mirror,valid)
      if(.not.valid) return
      reverse=[next_node,next_cell,edge(5),-edge(6)];key=edge
      if(route_edge_less(reverse,key)) key=reverse
    end subroutine canonical_route_edge
    pure logical function route_edge_less(left,right) result(less)
      integer, intent(in) :: left(6),right(6)
      integer :: component
      less=.false.
      do component=1,6
        if(left(component)==right(component)) cycle
        less=left(component)<right(component);return
      enddo
    end function route_edge_less
    recursive logical function regular_route_endpoint(node,cell) result(regular)
      integer, intent(in) :: node,cell(3)
      integer :: descendant,lo(3),hi(3),first(3),last(3),r,ii,jj,kk,a,face(3)
      regular=.false.
      if(geometry(node)%volume_fraction(cell(1),cell(2),cell(3))/=1) return
      do a=1,3
        face=cell;face(a)=face(a)-1
        if(open_area(node,a,face)/=open_area_full(node,a)) return
        if(open_area(node,a,cell)/=open_area_full(node,a)) return
      enddo
      descendant=covering_child(node,cell)
      if(descendant/=0) then
        call bounds(tree,descendant,lo,hi)
        r=tree%nodes(descendant)%patch%refinement_ratio
        first=(cell-lo)*r+1;last=first+r-1
        do kk=first(3),last(3)
        do jj=first(2),last(2)
        do ii=first(1),last(1)
          if(.not.regular_route_endpoint(descendant,[ii,jj,kk])) return
        enddo
        enddo
        enddo
      endif
      regular=.true.
    end function regular_route_endpoint
    real(dp) function open_area_full(node,a) result(area)
      integer, intent(in) :: node,a
      select case(a)
      case(1);area=geometry(node)%dy*geometry(node)%dz
      case(2);area=geometry(node)%dx*geometry(node)%dz
      case(3);area=geometry(node)%dx*geometry(node)%dy
      end select
    end function open_area_full
    ! Fully regular endpoint spreading is uniform in every direction.
    ! A failed endpoint cannot be repaired by changing its incident edge.
    subroutine exclude_regular_endpoint_reroute(node,cell,target,regular)
      integer, intent(in) :: node,cell(3),target
      logical, intent(in) :: regular(2)
      integer :: endpoint,root_node,root_cell(3),ancestor,index(3),lo(3),hi(3),parent_node,r
      do endpoint=1,2
        root_node=child;root_cell=positions(:,target)
        if(endpoint==2) then
          root_node=donor_node;root_cell=donor_cell
        endif
        if(.not.regular(endpoint)) cycle
        ancestor=node;index=cell
        do while(ancestor/=root_node)
          parent_node=tree%nodes(ancestor)%parent
          if(parent_node==0) exit
          call bounds(tree,ancestor,lo,hi)
          r=tree%nodes(ancestor)%patch%refinement_ratio
          index=lo+(index-1)/r;ancestor=parent_node
        enddo
        if(ancestor==root_node.and.all(index==root_cell)) route_can_reroute=.false.
      enddo
    end subroutine exclude_regular_endpoint_reroute
    subroutine recipient_allowance(target,maximum,valid,deficit)
      integer, intent(in) :: target
      real(dp), intent(out) :: maximum
      logical, intent(out) :: valid
      real(dp), optional, intent(out) :: deficit(nv)
      real(dp) :: low(nv),high(nv)
      low=0;high=0;maximum=1
      if(present(deficit)) deficit=0
      if(side==0) low=transfers(:,target)-weights(target)*old
      if(side==1) high=weights(target)*old-transfers(:,target)
      call propagate_cell_packet(child,positions(:,target),axis,low,high,.true.,valid,maximum,deficit)
    end subroutine recipient_allowance
    subroutine build_straight_route(target,length,a,step,valid)
      integer, intent(in) :: target,length,a,step
      logical, intent(out) :: valid
      integer :: entry,node,cell(3),next_node,next_cell(3),face(3),next_face(3),ios
      logical :: ready
      valid=.false.
      if(allocated(route_node)) deallocate(route_node,route_cell,route_axis,route_step)
      allocate(route_node(length),route_cell(3,length),route_axis(length),route_step(length),stat=ios)
      ready=ios==0
      if(.not.ready) return
      route_length=length;node=child;cell=positions(:,target)
      do entry=1,length
        route_node(entry)=node;route_cell(:,entry)=cell;route_axis(entry)=a;route_step(entry)=step
        call path_neighbor(node,cell,a,step,next_node,next_cell,face,next_face,valid)
        if(.not.valid) return
        node=next_node;cell=next_cell
      enddo
      valid=.true.
    end subroutine build_straight_route
    ! Test the complete donor/recipient/interior change before selecting a
    ! route. Trial descendants and faces are restored on either outcome, so
    ! a rejected route does not consume inventories needed by later searches.
    recursive subroutine route_accepts_packet(target,valid,allow_partial)
      integer, intent(in) :: target
      logical, intent(out) :: valid
      logical, optional, intent(in) :: allow_partial
      type(amr_tree_state_3d) :: saved_candidate
      type(eb_integrated_flux_3d), allocatable :: saved_working(:),saved_effective(:)
      logical, allocatable :: saved_changed(:)
      real(dp) :: primitive(reactive_nprim(size(species))),recovered,speed,value(nv),lower,upper,maximum, &
      baseline(nv),direction(nv)
      real(dp) :: before_allowance,after_allowance,prior_allowance,before_deficit(nv),after_deficit(nv)
      logical :: allowance_ok,endpoint_regular(2)
      integer :: node,cell(3),d(3),ii,jj,kk,entry,prior
      logical :: good,local_ok,needs_trial,propagation_ok,interval_ok
      route_complete=.false.;route_can_reroute=.false.
      route_improved=covering_child(child,positions(:,target))==0
      before_allowance=0
      if(.not.route_improved) then
        call recipient_allowance(target,before_allowance,allowance_ok,before_deficit)
        if(.not.allowance_ok) then
          valid=.false.;return
        endif
      endif
      local_ok=.true.
      do entry=1,2
        node=child;cell=positions(:,target);value=funded_fine(:,target)
        if(entry==2) then
          node=donor_node;cell=donor_cell;value=funded_donor(:,target)
        endif
        if(.not.all(ieee_is_finite(value)).or.value(1)<=0.or.any(value(6:active_components)<0)) then
          local_ok=.false.;cycle
        endif
        call recover_repair_cell(species,geometry(node),cell,value, &
          temperature%fields(node)%values(1,cell(1),cell(2),cell(3)),primitive,recovered,speed,good,frame, &
          repair_node=node)
        local_ok=local_ok.and.good
      enddo
      valid=local_ok
      if(.not.valid) then
        if(present(allow_partial)) then
          if(allow_partial) call find_partial_route_packet(target,1._dp,valid)
        endif
        return
      endif
      needs_trial=covering_child(donor_node,donor_cell)/=0.or.covering_child(child,positions(:,target))/=0
      do entry=2,route_length
        if(covering_child(route_node(entry),route_cell(:,entry))/=0) needs_trial=.true.
      enddo
      ! Only geometrically different endpoint spreading or a covered
      ! interior can change acceptance when the donor/path is rerouted.
      endpoint_regular(1)=regular_route_endpoint(child,positions(:,target))
      endpoint_regular(2)=regular_route_endpoint(donor_node,donor_cell)
      route_can_reroute=(covering_child(child,positions(:,target))/=0.and..not.endpoint_regular(1)).or. &
        (covering_child(donor_node,donor_cell)/=0.and..not.endpoint_regular(2))
      do entry=2,route_length
        if(covering_child(route_node(entry),route_cell(:,entry))/=0) route_can_reroute=.true.
      enddo
      ! Ordinary same-level routes leave every interior state unchanged.
      if(.not.needs_trial) return
      saved_candidate=candidate;saved_working=working;saved_effective=effective;saved_changed=changed
      cell=positions(:,target)
      if(covering_child(child,cell)==0) then
        candidate%fields(child)%values(:,cell(1),cell(2),cell(3))=funded_fine(:,target)
        changed(child)=.true.
      endif
      if(covering_child(donor_node,donor_cell)==0) then
        candidate%fields(donor_node)%values(:,donor_cell(1),donor_cell(2),donor_cell(3))=funded_donor(:,target)
        changed(donor_node)=.true.
      endif
      call propagate_route_interior(exchange(:,target),good)
      local_ok=good;propagation_ok=good;maximum=1
      if(good.and.covering_child(child,positions(:,target))/=0) then
        call recipient_allowance(target,after_allowance,allowance_ok,after_deficit)
        ! A different species can remain the limiting constraint while this
        ! packet reduces another deficit. Accept that complementary progress
        ! only when no component's missing physical inventory increases.
        route_improved=allowance_ok.and.(after_allowance>before_allowance.or. &
          (all(after_deficit<=before_deficit).and.any(after_deficit<before_deficit)))
      endif
      ! A route can redistribute descendants of an earlier recipient.
      ! Preserve its funded correction budget before spending this packet.
      if(good) then
        do prior=1,target-1
          call recipient_allowance(prior,prior_allowance,good)
          if(.not.good.or.prior_allowance<required_correction) local_ok=.false.
        enddo
      endif
      if(propagation_ok) then
        do node=1,size(tree%nodes)
          if(.not.changed(node)) cycle
          d=tree%nodes(node)%shape
          do kk=1,d(3)
          do jj=1,d(2)
          do ii=1,d(1)
            if(geometry(node)%volume_fraction(ii,jj,kk)<=0) cycle
            associate(trial_value=>candidate%fields(node)%values(:,ii,jj,kk))
              if(all(trial_value==saved_candidate%fields(node)%values(:,ii,jj,kk))) cycle
              baseline=saved_candidate%fields(node)%values(:,ii,jj,kk);direction=trial_value-baseline
              lower=0;upper=1
              call reactive_eb_exchange_interval_3d(reshape(baseline,[nv,1]),reshape(direction,[nv,1]), &
                lower,upper,interval_ok, &
        passive_components=passive_components)
              if(.not.interval_ok.or.lower>0) upper=0
              maximum=min(maximum,upper)
              if(upper<1) call exclude_regular_endpoint_reroute(node,[ii,jj,kk],target,endpoint_regular)
              if(.not.all(ieee_is_finite(trial_value)).or.trial_value(1)<=0.or.any(trial_value(6:active_components)<0)) then
                local_ok=.false.;cycle
              endif
              call recover_repair_cell(species,geometry(node),[ii,jj,kk],trial_value, &
                temperature%fields(node)%values(1,ii,jj,kk),primitive,recovered,speed,good,frame, &
          repair_node=node)
              local_ok=local_ok.and.good
              if(.not.good) call exclude_regular_endpoint_reroute(node,[ii,jj,kk],target,endpoint_regular)
            end associate
          enddo
          enddo
          enddo
        enddo
      endif
      good=local_ok
      if(good.and.covering_child(child,positions(:,target))/=0) then
        if(allowance_ok.and.after_allowance>=required_correction) then
          do prior=1,target
            call project_recipient_endpoint(prior,good)
            if(.not.good) exit
          enddo
          if(good) call accept_route_candidate(saved_candidate,good)
          local_ok=good;route_complete=good
        endif
      endif
      ! Propagation failures may leave a partially populated private ledger.
      if(allocated(packet_keys)) deallocate(packet_keys)
      if(allocated(packet_values)) deallocate(packet_values)
      candidate=saved_candidate;working=saved_working;effective=saved_effective;changed=saved_changed
      valid=local_ok
      ! A parent average may admit more withdrawal than its least stocked
      ! descendant. In the pooling search, retain the physical fraction and
      ! verify its represented states and EOS again before spending it.
      if(valid.or..not.propagation_ok.or.maximum<=0) return
      if(.not.present(allow_partial)) return
      if(.not.allow_partial) return
      call find_partial_route_packet(target,maximum,valid)
    end subroutine route_accepts_packet
    ! Retain donor inventory once the covered recipient can accept its final
    ! interface correction. Each trial uses the same complete conservative
    ! direction and validates physical descendants plus the corrected EOS.
    subroutine minimize_recipient_packet(target,valid)
      integer, intent(in) :: target
      logical, intent(out) :: valid
      real(dp) :: original_fine(nv),original_donor(nv),original_exchange(nv)
      real(dp) :: lower,upper,middle
      integer :: iteration
      logical :: accepted
      original_fine=funded_fine(:,target);original_donor=funded_donor(:,target);original_exchange=exchange(:,target)
      lower=0;upper=1
      do iteration=1,64
        middle=lower+(upper-lower)/2
        if(middle<=lower.or.middle>=upper) exit
        exchange(:,target)=middle*original_exchange
        funded_fine(:,target)=fine_state(:,target)+exchange(:,target)/volumes(target)
        funded_donor(:,target)=donor_state(:,target)-exchange(:,target)/donor_volume(target)
        call route_accepts_packet(target,accepted)
        if(accepted.and.route_complete) then
          upper=middle
        else
          lower=middle
        endif
      enddo
      ! Retain a represented margin for the later group reconciliation,
      ! whose component sums may round differently from this endpoint trial.
      upper=min(1._dp,upper+128*epsilon(1._dp))
      exchange(:,target)=upper*original_exchange
      funded_fine(:,target)=fine_state(:,target)+exchange(:,target)/volumes(target)
      funded_donor(:,target)=donor_state(:,target)-exchange(:,target)/donor_volume(target)
      call route_accepts_packet(target,valid)
      valid=valid.and.route_complete.and.route_improved
      if(valid) return
      ! Preserve the already accepted full packet if represented endpoint
      ! rounding prevents acceptance of the reduced packet.
      funded_fine(:,target)=original_fine;funded_donor(:,target)=original_donor;exchange(:,target)=original_exchange
      call route_accepts_packet(target,valid)
      valid=valid.and.route_complete.and.route_improved
    end subroutine minimize_recipient_packet
    subroutine project_recipient_endpoint(target,valid)
      integer, intent(in) :: target
      logical, intent(out) :: valid
      real(dp) :: low(nv),high(nv)
      integer :: count,ios
      logical :: ready
      valid=.false.;count=transfer_cell_count(child,positions(:,target))
      if(count<0.or.count>huge(count)/2) return
      allocate(packet_keys(5,2*count),packet_values(nv,2*count),stat=ios)
      ready=ios==0
      if(.not.ready) return
      packet_count=0;low=0;high=0
      if(side==0) low=required_correction*(transfers(:,target)-weights(target)*old)
      if(side==1) high=required_correction*(weights(target)*old-transfers(:,target))
      call propagate_cell_packet(child,positions(:,target),axis,low,high,.true.,valid)
      deallocate(packet_keys,packet_values)
    end subroutine project_recipient_endpoint
    subroutine accept_route_candidate(baseline,valid)
      type(amr_tree_state_3d), intent(in) :: baseline
      logical, intent(out) :: valid
      integer :: node,d(3),ii,jj,kk
      real(dp) :: primitive(reactive_nprim(size(species))),recovered,speed
      logical :: good,local_ok
      local_ok=.true.
      do node=1,size(tree%nodes)
        if(.not.changed(node)) cycle
        d=tree%nodes(node)%shape
        do kk=1,d(3)
        do jj=1,d(2)
        do ii=1,d(1)
          if(geometry(node)%volume_fraction(ii,jj,kk)<=0) cycle
          associate(value=>candidate%fields(node)%values(:,ii,jj,kk))
            if(all(value==baseline%fields(node)%values(:,ii,jj,kk))) cycle
            if(.not.all(ieee_is_finite(value)).or.value(1)<=0.or.any(value(6:active_components)<0)) then
              local_ok=.false.;cycle
            endif
            call recover_repair_cell(species,geometry(node),[ii,jj,kk],value, &
              temperature%fields(node)%values(1,ii,jj,kk),primitive,recovered,speed,good,frame, &
          repair_node=node)
            local_ok=local_ok.and.good
          end associate
        enddo
        enddo
        enddo
      enddo
      valid=local_ok
    end subroutine accept_route_candidate
    ! Keep all conserved components on the same exchange direction. The
    ! density/species bound supplies an upper trial; EOS recovery then brackets
    ! a represented admissible fraction without altering its tolerances.
    recursive subroutine find_partial_route_packet(target,maximum,valid)
      integer, intent(in) :: target
      real(dp), intent(in) :: maximum
      logical, intent(out) :: valid
      real(dp) :: original_fine(nv),original_donor(nv),original_exchange(nv)
      real(dp) :: lower,upper,middle
      integer :: iteration
      logical :: accepted
      valid=.false.
      original_fine=funded_fine(:,target);original_donor=funded_donor(:,target);original_exchange=exchange(:,target)
      upper=maximum*(1-32*epsilon(1._dp));lower=0
      exchange(:,target)=(upper)*original_exchange
      funded_fine(:,target)=fine_state(:,target)+exchange(:,target)/volumes(target)
      funded_donor(:,target)=donor_state(:,target)-exchange(:,target)/donor_volume(target)
      call route_accepts_packet(target,accepted)
      if(accepted.and.any(exchange(:,target)/=0)) then
        valid=.true.;return
      endif
      exchange(:,target)=(0._dp)*original_exchange
      funded_fine(:,target)=fine_state(:,target)+exchange(:,target)/volumes(target)
      funded_donor(:,target)=donor_state(:,target)-exchange(:,target)/donor_volume(target)
      call route_accepts_packet(target,accepted)
      if(accepted) then
        do iteration=1,64
          middle=lower+(upper-lower)/2
          if(middle<=lower.or.middle>=upper) exit
          exchange(:,target)=(middle)*original_exchange
          funded_fine(:,target)=fine_state(:,target)+exchange(:,target)/volumes(target)
          funded_donor(:,target)=donor_state(:,target)-exchange(:,target)/donor_volume(target)
          call route_accepts_packet(target,accepted)
          if(accepted) then
            lower=middle
          else
            upper=middle
          endif
        enddo
        if(lower>0) then
          exchange(:,target)=(lower)*original_exchange
          funded_fine(:,target)=fine_state(:,target)+exchange(:,target)/volumes(target)
          funded_donor(:,target)=donor_state(:,target)-exchange(:,target)/donor_volume(target)
          call route_accepts_packet(target,valid)
          valid=valid.and.any(exchange(:,target)/=0)
        endif
      endif
      if(.not.valid) then
        funded_fine(:,target)=original_fine;funded_donor(:,target)=original_donor;exchange(:,target)=original_exchange
      endif
    end subroutine find_partial_route_packet
    ! A turning route can change descendant inventories even when its
    ! parent cell has zero net divergence. Propagate each axis once, with both
    ! incident faces assembled, and retain the transaction's final EOS gate.
    subroutine propagate_route_interior(packet,valid)
      real(dp), intent(in) :: packet(nv)
      logical, intent(out) :: valid
      integer :: entry,node,cell(3),a,count,capacity,ios
      real(dp) :: low(nv),high(nv)
      logical :: ready
      valid=.false.;capacity=0
      do entry=2,route_length
        node=route_node(entry);cell=route_cell(:,entry)
        if(covering_child(node,cell)==0) cycle
        count=transfer_cell_count(node,cell)
        if(count<0.or.count>(huge(capacity)-capacity)/6) return
        capacity=capacity+6*count
      enddo
      if(covering_child(child,route_cell(:,1))/=0) then
        count=transfer_cell_count(child,route_cell(:,1))
        if(count<0.or.count>(huge(capacity)-capacity)/2) return
        capacity=capacity+2*count
      endif
      if(covering_child(donor_node,donor_cell)/=0) then
        count=transfer_cell_count(donor_node,donor_cell)
        if(count<0.or.count>(huge(capacity)-capacity)/2) return
        capacity=capacity+2*count
      endif
      if(capacity==0) then
        valid=.true.;return
      endif
      allocate(packet_keys(5,capacity),packet_values(nv,capacity),stat=ios)
      ready=ios==0
      if(.not.ready) return
      packet_count=0
      if(covering_child(child,route_cell(:,1))/=0) then
        a=route_axis(1);low=0;high=0
        if(route_step(1)>0) high=-packet
        if(route_step(1)<0) low=packet
        call propagate_cell_packet(child,route_cell(:,1),a,low,high,.true.,valid)
        if(.not.valid) return
      endif
      do entry=2,route_length
        node=route_node(entry);cell=route_cell(:,entry)
        if(covering_child(node,cell)==0) cycle
        do a=1,3
          low=0;high=0
          if(route_axis(entry-1)==a) then
            if(route_step(entry-1)>0) low=-packet
            if(route_step(entry-1)<0) high=packet
          endif
          if(route_axis(entry)==a) then
            if(route_step(entry)>0) high=-packet
            if(route_step(entry)<0) low=packet
          endif
          if(all(low==0).and.all(high==0)) cycle
          call propagate_cell_packet(node,cell,a,low,high,.true.,valid)
          if(.not.valid) return
        enddo
      enddo
      if(covering_child(donor_node,donor_cell)/=0) then
        ! The last edge enters the donor. Its incident packet withdraws the
        ! exchange from every physical descendant and updates their faces.
        a=route_axis(route_length);low=0;high=0
        if(route_step(route_length)>0) low=-packet
        if(route_step(route_length)<0) high=packet
        call propagate_cell_packet(donor_node,donor_cell,a,low,high,.true.,valid)
        if(.not.valid) return
      endif
      deallocate(packet_keys,packet_values)
      valid=.true.
    end subroutine propagate_route_interior
    integer function covering_child(node,cell) result(covered_by)
      integer, intent(in) :: node,cell(3)
      integer :: peer,lower(3),upper(3)
      covered_by=0
      do peer=1,size(tree%nodes)
        if(tree%nodes(peer)%parent/=node) cycle
        call bounds(tree,peer,lower,upper)
        if(all(cell>=lower).and.all(cell<=upper)) then
          covered_by=peer;return
        endif
      enddo
    end function covering_child
    recursive integer function transfer_cell_count(node,cell) result(count)
      integer, intent(in) :: node,cell(3)
      integer :: covered_by,lower(3),upper(3),first(3),last(3),r,ii,jj,kk,addition
      count=1;covered_by=covering_child(node,cell)
      if(covered_by==0) return
      call bounds(tree,covered_by,lower,upper)
      r=tree%nodes(covered_by)%patch%refinement_ratio
      first=(cell-lower)*r+1;last=first+r-1
      do kk=first(3),last(3)
      do jj=first(2),last(2)
      do ii=first(1),last(1)
        if(cell_volume(covered_by,[ii,jj,kk])<=0) cycle
        addition=transfer_cell_count(covered_by,[ii,jj,kk])
        if(addition<0.or.count>huge(count)-addition) then
          count=-1;return
        endif
        count=count+addition
      enddo
      enddo
      enddo
    end function transfer_cell_count
    real(dp) function cell_volume(node,cell) result(volume)
      integer, intent(in) :: node,cell(3)
      volume=geometry(node)%volume_fraction(cell(1),cell(2),cell(3))* &
        geometry(node)%dx*geometry(node)%dy*geometry(node)%dz
    end function cell_volume
    ! Record each shared internal face once. Adjacent refined cells supply the
    ! same packet to their common face; a disagreement aborts the transaction.
    subroutine record_refined_packet(node,a,face,packet,valid)
      integer, intent(in) :: node,a,face(3)
      real(dp), intent(in) :: packet(nv)
      logical, intent(out) :: valid
      integer :: entry,key(5)
      valid=.true.
      if(all(packet==0)) return
      if(open_area(node,a,face)<=0) then
        valid=.false.;return
      endif
      key=[node,a,face]
      do entry=1,packet_count
        if(any(packet_keys(:,entry)/=key)) cycle
        valid=all(packet_values(:,entry)==packet);return
      enddo
      if(packet_count==size(packet_keys,2)) then
        valid=.false.;return
      endif
      packet_count=packet_count+1;packet_keys(:,packet_count)=key;packet_values(:,packet_count)=packet
      call add_path_packet(node,a,face,packet)
    end subroutine record_refined_packet
    recursive subroutine propagate_cell_packet(node,cell,a,low_packet,high_packet,root,valid,maximum_correction,deficit)
      integer, intent(in) :: node,cell(3),a
      real(dp), intent(in) :: low_packet(nv),high_packet(nv)
      logical, intent(in) :: root
      logical, intent(out) :: valid
      real(dp), optional, intent(inout) :: maximum_correction,deficit(nv)
      integer :: covered_by,lower(3),upper(3),first(3),last(3),index(3),face(3),r,ii,jj,kk,along
      real(dp) :: volume,total_volume,column_volume,used_volume,low_area,high_area,all_low_area,all_high_area
      real(dp) :: incoming(nv),outgoing(nv),previous(nv),next(nv),budget(nv),slope(nv),delta(nv),exchange_lower,exchange_upper
      valid=.false.;volume=cell_volume(node,cell)
      if(volume<=0) return
      if(present(maximum_correction)) then
        budget=candidate%fields(node)%values(:,cell(1),cell(2),cell(3))
        if(.not.all(ieee_is_finite(budget)).or.budget(1)<=0.or.any(budget(6:active_components)<0)) return
        slope=(low_packet-high_packet)/volume
        if(.not.all(ieee_is_finite(slope))) return
        exchange_lower=0;exchange_upper=1
        call reactive_eb_exchange_interval_3d(reshape(budget,[nv,1]),reshape(slope,[nv,1]), &
          exchange_lower,exchange_upper,valid, &
        passive_components=passive_components)
        if(.not.valid.or.exchange_lower>0) return
        maximum_correction=min(maximum_correction,exchange_upper)
        if(tree%nodes(parent)%level<minimum_level.and.tree%nodes(node)%level<=minimum_level) then
          call cell_increment(node,cell,[0,0,0],0,delta)
          call add_eb_tree_delta_cell_3d(budget,delta,valid)
          if(.not.valid) return
          exchange_lower=0;exchange_upper=1
          call reactive_eb_exchange_interval_3d(reshape(budget,[nv,1]),reshape(slope,[nv,1]), &
            exchange_lower,exchange_upper,valid,passive_components=passive_components)
          if(.not.valid.or.exchange_lower>0) return
          maximum_correction=min(maximum_correction,exchange_upper)
        endif
        if(.not.root) then
          face=cell;face(a)=face(a)-1
          if(open_area(node,a,face)<=0.and.any(low_packet/=0)) return
          if(open_area(node,a,cell)<=0.and.any(high_packet/=0)) return
        endif
      else
        candidate%fields(node)%values(:,cell(1),cell(2),cell(3))= &
          candidate%fields(node)%values(:,cell(1),cell(2),cell(3))+(low_packet-high_packet)/volume
        changed(node)=.true.
        if(.not.root) then
          face=cell;face(a)=face(a)-1
          call record_refined_packet(node,a,face,low_packet,valid)
          if(.not.valid) return
          call record_refined_packet(node,a,cell,high_packet,valid)
          if(.not.valid) return
        endif
      endif
      covered_by=covering_child(node,cell)
      valid=.true.
      if(covered_by==0) then
        ! Count physical leaves once; covered averages are not inventories.
        ! MPI broadcasts above keep this traversal identical on every rank.
        if(present(deficit)) then
          if(.not.present(maximum_correction)) then
            valid=.false.;return
          endif
          budget=budget+required_correction*slope
          if(.not.all(ieee_is_finite(budget))) then
            valid=.false.;return
          endif
          deficit(1)=deficit(1)+volume*max(0._dp,-budget(1))
          deficit(6:active_components)=deficit(6:active_components)+volume*max(0._dp,-budget(6:active_components))
          valid=all(ieee_is_finite(deficit))
        endif
        return
      endif
      call bounds(tree,covered_by,lower,upper)
      r=tree%nodes(covered_by)%patch%refinement_ratio
      first=(cell-lower)*r+1;last=first+r-1
      total_volume=0;all_low_area=0;all_high_area=0
      do kk=first(3),last(3)
      do jj=first(2),last(2)
      do ii=first(1),last(1)
        index=[ii,jj,kk];total_volume=total_volume+cell_volume(covered_by,index)
        if(index(a)/=first(a)) cycle
        face=index;face(a)=first(a)-1;all_low_area=all_low_area+open_area(covered_by,a,face)
        face(a)=last(a);all_high_area=all_high_area+open_area(covered_by,a,face)
      enddo
      enddo
      enddo
      valid=abs(total_volume-volume)<=64*epsilon(1._dp)*volume
      if(.not.valid) return
      if((all_low_area<=0.and.any(low_packet/=0)).or.(all_high_area<=0.and.any(high_packet/=0))) then
        valid=.false.;return
      endif
      do kk=first(3),last(3)
      do jj=first(2),last(2)
      do ii=first(1),last(1)
        index=[ii,jj,kk]
        if(index(a)/=first(a)) cycle
        face=index;face(a)=first(a)-1;low_area=open_area(covered_by,a,face)
        face(a)=last(a);high_area=open_area(covered_by,a,face)
        incoming=0;outgoing=0
        if(all_low_area>0) incoming=(low_area/all_low_area)*low_packet
        if(all_high_area>0) outgoing=(high_area/all_high_area)*high_packet
        column_volume=0
        do along=first(a),last(a)
          index(a)=along;column_volume=column_volume+cell_volume(covered_by,index)
        enddo
        if(column_volume<=0) then
          valid=all(incoming==0).and.all(outgoing==0)
          if(.not.valid) return
          cycle
        endif
        used_volume=0;previous=incoming
        do along=first(a),last(a)
          index(a)=along;used_volume=used_volume+cell_volume(covered_by,index)
          if(along==last(a)) then
            next=outgoing
          else
            next=(1-used_volume/column_volume)*incoming+(used_volume/column_volume)*outgoing
          endif
          if(cell_volume(covered_by,index)>0) then
            call propagate_cell_packet(covered_by,index,a,previous,next,.false.,valid,maximum_correction,deficit)
            if(.not.valid) return
          else if(any(previous/=next)) then
            valid=.false.;return
          endif
          previous=next
        enddo
      enddo
      enddo
      enddo
    end subroutine propagate_cell_packet
    ! Locate a straight path in the integer layout, including periodic seams.
    ! Coarse gaps and solid cells stop a route. Covered cells carry or donate
    ! through their physical descendants, with the same trial acceptance gate.
    ! Both copies of a patch-boundary face carry the same extensive packet.
    subroutine path_neighbor(node,cell,a,step,next_node,next_cell,face,next_face,valid)
      integer, intent(in) :: node,cell(3),a,step
      integer, intent(out) :: next_node,next_cell(3),face(3),next_face(3)
      logical, intent(out) :: valid
      integer(int64) :: global_cell(3)
      integer :: peer,depth
      valid=.false.;next_node=0;next_cell=0;face=cell;face(a)=face(a)+min(step,0)
      depth=tree%nodes(node)%level+1
      if(.not.active_path_cell(node,cell)) return
      global_cell=origin(:,node)+int(cell-1,int64);global_cell(a)=global_cell(a)+step
      if(periodic(a)) global_cell(a)=modulo(global_cell(a),domain(a,depth))
      if(any(global_cell<0).or.any(global_cell>=domain(:,depth))) return
      do peer=1,size(tree%nodes)
        if(tree%nodes(peer)%level/=depth-1) cycle
        if(any(global_cell<origin(:,peer)).or.any(global_cell>=origin(:,peer)+tree%nodes(peer)%shape)) cycle
        next_node=peer;next_cell=int(global_cell-origin(:,peer))+1;exit
      enddo
      if(next_node==0) return
      if(.not.active_path_cell(next_node,next_cell)) return
      next_face=next_cell;next_face(a)=next_face(a)+min(-step,0)
      valid=open_area(node,a,face)>0.and.open_area(next_node,a,next_face)>0
    end subroutine path_neighbor
    logical function active_path_cell(node,cell) result(active)
      integer, intent(in) :: node,cell(3)
      active=geometry(node)%volume_fraction(cell(1),cell(2),cell(3))>0
    end function active_path_cell
    subroutine add_path_packet(node,a,face,packet)
      integer, intent(in) :: node,a,face(3)
      real(dp), intent(in) :: packet(nv)
      real(dp) :: intensive(nv)
      intensive=packet/open_area(node,a,face)
      call set_face(working(node),a,face,face_value(working(node),a,face)+intensive)
      call set_face(effective(node),a,face,face_value(effective(node),a,face)+intensive)
    end subroutine add_path_packet
    real(dp) function open_area(p,a,index) result(value)
      integer, intent(in) :: p,a,index(3)
      associate(g=>geometry(p))
      select case(a)
      case(1);value=g%dy*g%dz*g%x_face_fraction(index(1),index(2),index(3))
      case(2);value=g%dx*g%dz*g%y_face_fraction(index(1),index(2),index(3))
      case(3);value=g%dx*g%dy*g%z_face_fraction(index(1),index(2),index(3))
      end select
      end associate
    end function open_area
    function face_value(f,a,index) result(value)
      type(eb_integrated_flux_3d), intent(in) :: f
      integer, intent(in) :: a,index(3)
      real(dp) :: value(nv)
      select case(a)
      case(1);value=f%x(:,index(1),index(2),index(3))
      case(2);value=f%y(:,index(1),index(2),index(3))
      case(3);value=f%z(:,index(1),index(2),index(3))
      end select
    end function face_value
    subroutine set_face(f,a,index,value)
      type(eb_integrated_flux_3d), intent(inout) :: f
      integer, intent(in) :: a,index(3)
      real(dp), intent(in) :: value(nv)
      select case(a)
      case(1);f%x(:,index(1),index(2),index(3))=value
      case(2);f%y(:,index(1),index(2),index(3))=value
      case(3);f%z(:,index(1),index(2),index(3))=value
      end select
    end subroutine set_face
    subroutine cell_increment(p,cell,excluded,excluded_axis,value)
      integer, intent(in) :: p,cell(3),excluded(3),excluded_axis
      real(dp), intent(out) :: value(nv)
      integer :: a,step,index(3)
      real(dp) :: volume,scale
      value=0
      if(present(rereflux_delta).and.excluded_axis==0) &
        value=routed_register%fields(p)%values(:,cell(1),cell(2),cell(3))
      if(immediate_funding) return
      volume=geometry(p)%volume_fraction(cell(1),cell(2),cell(3))*geometry(p)%dx*geometry(p)%dy*geometry(p)%dz
      do a=1,3
      do step=-1,1,2
        index=cell;index(a)=index(a)+min(step,0)
        if(a==excluded_axis.and.all(index==excluded)) cycle
        scale=-real(step,dp)*open_area(p,a,index)/volume
        value=value+scale*(face_value(effective(p),a,index)-face_value(working(p),a,index))
      enddo
      enddo
    end subroutine cell_increment
  end subroutine repair_reactive_eb_tree_interface_components_3d

  ! Input states have already accepted each node's raw time-integrated flux.
  ! Replace shared transfers from finest to coarsest on exposed cells only,
  ! then average down. Publish corrected states and effective fluxes together.
  ! Close one accepted interface interval. Only this interface's coarse
  ! and fine ledgers are consumed, after gas/flux EOS acceptance succeeds.
  subroutine synchronize_reactive_eb_tree_register_flux_3d(species,tree,flux,extended_geometry,registers,depth,ok, &
      boundary_faces,frame)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_reactive_tree_3d), intent(inout) :: tree
    type(eb_integrated_flux_3d), intent(inout) :: flux(:)
    type(eb_geometry_3d), intent(in) :: extended_geometry(:)
    type(state_redist_register_3d), intent(inout) :: registers(:)
    integer, intent(in) :: depth
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d) :: correction,spill
    type(eb_reactive_tree_3d) :: candidate
    type(eb_integrated_flux_3d), allocatable :: working(:)
    integer :: p,nv,first_depth
    logical :: good,has_spill
    ok=.false.;nv=reactive_nvar(size(species))
    if(.not.valid_tree(tree,nv)) return
    call assemble_state_redist_interface_rereflux_3d(tree%topology,extended_geometry,registers,nv,depth, &
      correction,good,boundary_faces,spill,has_spill)
    if(.not.good) return
    first_depth=depth
    if(has_spill) then
      first_depth=depth-1
      do p=1,size(tree%nodes)
        if(tree%topology%nodes(p)%level==first_depth) correction%fields(p)%values=spill%fields(p)%values
      enddo
    endif
    candidate=tree;working=flux
    call synchronize_reactive_eb_tree_flux_3d(species,candidate,working,good,depth,depth+1, &
      boundary_faces,frame,correction,rereflux_first_level=first_depth)
    if(.not.good) then
      candidate=tree;working=flux
      call repair_reactive_eb_tree_interface_budget_3d(species,candidate,working,good,depth, &
        boundary_faces,frame,correction,rereflux_first_level=first_depth)
      if(.not.good) return
      call synchronize_reactive_eb_tree_flux_3d(species,candidate,working,good,depth,depth+1, &
        boundary_faces,frame,correction,rereflux_first_level=first_depth)
      if(.not.good) return
    end if
    tree=candidate;flux=working
    do p=1,size(tree%nodes)
      if(tree%topology%nodes(p)%level==depth) then
        registers(p)%coarse_delta=0
      else if(tree%topology%nodes(p)%level==depth+1) then
        registers(p)%fine_delta=0
      endif
    enddo
    ok=.true.
  end subroutine synchronize_reactive_eb_tree_register_flux_3d

  subroutine synchronize_reactive_eb_tree_flux_3d(species,tree,flux,ok,minimum_level,maximum_level, &
      boundary_faces,frame,rereflux_delta,rereflux_first_level)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_reactive_tree_3d), intent(inout) :: tree
    type(eb_integrated_flux_3d), intent(inout) :: flux(:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: minimum_level,maximum_level
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d), optional, intent(in) :: rereflux_delta
    integer, optional, intent(in) :: rereflux_first_level
    call synchronize_reactive_eb_tree_flux_core_3d(species,tree,flux,ok,minimum_level,maximum_level, &
      boundary_faces,frame,rereflux_delta=rereflux_delta,rereflux_first_level=rereflux_first_level)
  end subroutine synchronize_reactive_eb_tree_flux_3d

  ! Signed composite reflux assembly, before EOS acceptance or averaging.
  ! The input tree and faces stay read-only; output publication is transactional.
  ! Readonly ancestor face debt. Pack gas and actual H into the existing
  ! conservative assembly so both channels use identical peer/CF ownership.
  ! Virtual coarse corrections remain available to the assembler, but never
  ! become physical inventory in the returned inactive donor stock.
  subroutine assemble_reactive_eb_held_funding_debt_3d(species,tree,stored,flux,scalar_flux,clock,debt,ok, &
      boundary_faces,frame,rereflux_delta,scalar_rereflux_delta)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: stored
    type(eb_integrated_flux_3d), intent(in) :: flux(:),scalar_flux(:)
    integer, intent(in) :: clock
    type(amr_tree_state_3d), intent(inout) :: debt
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(eb_reactive_tree_3d) :: packed_tree
    type(eb_geometry_3d), allocatable :: geometry(:)
    type(amr_tree_state_3d), optional, intent(in) :: rereflux_delta,scalar_rereflux_delta
    type(amr_tree_state_3d), allocatable :: joint_rereflux
    type(amr_tree_state_3d) :: raw,physical
    type(eb_integrated_flux_3d), allocatable :: packed_faces(:),canonical(:)
    logical, allocatable :: local(:)
    integer :: nv,p,n(3),status
    logical :: good
    ok=.false.;nv=reactive_nvar(size(species))
    if(.not.valid_tree(tree,nv)) return
    if(clock<0.or.clock>=tree%topology%level_count()) return
    if(.not.valid_held_mol_field_3d(tree,stored,1)) return
    allocate(geometry(size(tree%nodes)),local(size(tree%nodes)),stat=status)
    if(status/=0) return
    local=.true.
    do p=1,size(tree%nodes)
      geometry(p)=tree%nodes(p)%geometry
    enddo
    if(present(rereflux_delta).neqv.present(scalar_rereflux_delta)) return
    if(present(rereflux_delta)) then
      allocate(joint_rereflux,stat=status);if(status/=0) return
      call pack_eb_held_funding_increment_3d(tree%topology,geometry,nv,local,rereflux_delta, &
        scalar_rereflux_delta,joint_rereflux,good)
      if(.not.good) return
    endif
    call pack_eb_held_funding_faces_3d(tree%topology,geometry,nv,local,flux,scalar_flux,packed_faces,good)
    if(.not.good) return
    packed_tree=tree
    do p=1,size(tree%nodes)
      n=tree%topology%nodes(p)%shape
      deallocate(packed_tree%nodes(p)%state)
      allocate(packed_tree%nodes(p)%state(nv+1,n(1),n(2),n(3)),stat=status)
      if(status/=0) return
      packed_tree%nodes(p)%state(:nv,:,:,:)=tree%nodes(p)%state
      packed_tree%nodes(p)%state(nv+1:nv+1,:,:,:)=stored%fields(p)%values
    enddo
    call assemble_reactive_eb_tree_reflux_3d(species,packed_tree,packed_faces,raw,canonical,good, &
      0,clock,boundary_faces,frame,rereflux_delta=joint_rereflux,passive_components=1)
    if(.not.good) return
    call physical_eb_held_funding_debt_3d(tree%topology,geometry,nv,local,clock,raw,physical,good)
    if(.not.good) return
    call move_alloc(physical%fields,debt%fields);ok=.true.
  end subroutine assemble_reactive_eb_held_funding_debt_3d

  subroutine assemble_reactive_eb_tree_reflux_3d(species,tree,flux,increment,replacement,ok, &
      minimum_level,maximum_level,boundary_faces,frame,rereflux_delta,rereflux_first_level, &
      passive_components)
    integer, optional, intent(in) :: passive_components
    type(nasa7_species), intent(in) :: species(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(eb_integrated_flux_3d), intent(in) :: flux(:)
    type(amr_tree_state_3d), intent(inout) :: increment
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: replacement(:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: minimum_level,maximum_level,rereflux_first_level
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d), optional, intent(in) :: rereflux_delta
    type(eb_reactive_tree_3d) :: private_tree
    type(amr_tree_state_3d) :: private_increment
    type(eb_integrated_flux_3d), allocatable :: private_flux(:)
    private_tree=tree;private_flux=flux
    call synchronize_reactive_eb_tree_flux_core_3d(species,private_tree,private_flux,ok, &
      minimum_level,maximum_level,boundary_faces,frame,private_increment,rereflux_delta,rereflux_first_level, &
        passive_components=passive_components)
    if(.not.ok) return
    call move_alloc(private_increment%fields,increment%fields)
    call move_alloc(private_flux,replacement)
  end subroutine assemble_reactive_eb_tree_reflux_3d

  subroutine synchronize_reactive_eb_tree_flux_core_3d(species,tree,flux,ok, &
      minimum_level,maximum_level,boundary_faces,frame,reflux_increment,rereflux_delta,rereflux_first_level, &
      passive_components)
    integer, optional, intent(in) :: passive_components
    integer :: passives
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(eb_reactive_tree_3d), intent(inout) :: tree
    type(eb_integrated_flux_3d), intent(inout) :: flux(:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: minimum_level,maximum_level
    character(len=*), optional, intent(in) :: boundary_faces(6)
    character(len=16) :: faces(6)
    real(dp) :: period(3),length(3)
    type(amr_tree_state_3d), optional, intent(inout) :: reflux_increment
    ! Already assembled coarse re-reflux density corrections, including
    ! fine-covered cells. The register/neighbor assembly belongs to the caller.
    type(amr_tree_state_3d), optional, intent(in) :: rereflux_delta
    integer, optional, intent(in) :: rereflux_first_level
    type(amr_tree_state_3d) :: candidate_increment,routed_delta
    type(eb_reactive_tree_3d) :: candidate
    type(eb_integrated_flux_3d), allocatable :: effective(:)
    type(eb_integrated_flux_3d) :: restricted,owned
    type(eb_geometry_3d) :: section
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer(int64) :: offset(3),section_origin(3)
    integer :: nv,node,peer,owner,receiver,axis,depth,parent,child,lo(3),hi(3),d(3),ratio,first_level,last_level, &
      register_first_level,seam,side,i,j,k
    logical, allocatable :: refined(:,:,:)
    logical :: valid
    ok=.false.;passives=0
    if(present(passive_components)) passives=passive_components
    if(passives<0.or.passives>1) return
    if(passives>0.and..not.present(reflux_increment)) return
    nv=reactive_nvar(size(species))+passives
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    if (.not.valid_tree(tree,nv)) return
    faces='outflow'
    if(present(boundary_faces)) faces=boundary_faces
    if(any(faces/='outflow'.and.faces/='inflow'.and.faces/='periodic'.and. &
      faces/='slip_wall'.and.faces/='no_slip_wall'.and.faces/='isothermal_wall')) return
    if(.not.periodic_eb_faces_match(tree%nodes(1)%geometry,faces)) return
    associate(g=>tree%nodes(1)%geometry)
    length=[g%x_upper-g%x_lower,g%y_upper-g%y_lower,g%z_upper-g%z_lower]
    end associate
    first_level=0;last_level=tree%topology%level_count()-1
    if (present(minimum_level)) first_level=minimum_level
    if (present(maximum_level)) last_level=maximum_level
    if (first_level<0.or.last_level<first_level.or.last_level>=tree%topology%level_count()) return
    register_first_level=first_level
    if(present(rereflux_first_level)) register_first_level=rereflux_first_level
    if(register_first_level<0.or.register_first_level>first_level) return
    if (size(flux)/=size(tree%nodes)) return
    do node=1,size(flux)
      if (.not.valid_eb_integrated_flux_3d(flux(node),tree%nodes(node)%geometry,nv)) return
    end do
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if (.not.valid) return
    candidate=tree;effective=flux
    if(present(rereflux_delta)) then
      call route_eb_tree_rereflux_3d(tree,rereflux_delta,routed_delta,nv,register_first_level,last_level,valid)
      if(.not.valid) return
    endif
    if(present(reflux_increment).or.present(rereflux_delta)) then
      allocate(candidate_increment%fields(size(tree%nodes)))
      do node=1,size(tree%nodes)
        d=tree%topology%nodes(node)%shape
        allocate(candidate_increment%fields(node)%values(nv,d(1),d(2),d(3)))
        candidate_increment%fields(node)%values=0
      enddo
    endif
    ! Reconcile raw peer faces first at every depth. Descendant replacements
    ! below are mirrored into both incident parent patches and must not be
    ! overwritten later by a coarser owner decision.
    do node=1,size(flux)
    do peer=node+1,size(flux)
      if (tree%topology%nodes(node)%level<first_level.or.tree%topology%nodes(node)%level>last_level) cycle
      if (tree%topology%nodes(peer)%level/=tree%topology%nodes(node)%level) cycle
      owner=node;receiver=peer
      do axis=1,3
        if (origin(axis,node)==origin(axis,peer)) cycle
        if (origin(axis,node)>origin(axis,peer)) then
          owner=peer;receiver=node
        end if
        exit
      end do
      offset=origin(:,receiver)-origin(:,owner)
      if (any(offset < -int(tree%topology%nodes(receiver)%shape,int64)).or. &
        any(offset>int(tree%topology%nodes(owner)%shape,int64))) cycle
      call reconcile_reactive_eb_peer_flux_3d(species,tree%nodes(owner)%geometry,tree%nodes(receiver)%geometry, &
        int(offset),effective(owner),effective(receiver),ok=valid, &
        component_count=nv)
      if (.not.valid) return
    end do
    end do
    ! Canonical raw periodic transfer comes from the high-domain face.
    ! Copy owner flux so a patch spanning the domain can receive its own seam.
    do seam=1,3
      if(faces(2*seam-1)/='periodic') cycle
      period=0;period(seam)=length(seam)
      do owner=1,size(flux)
        depth=tree%topology%nodes(owner)%level
        if(depth<first_level.or.depth>last_level) cycle
        if(origin(seam,owner)+tree%topology%nodes(owner)%shape(seam)/=domain(seam,depth+1)) cycle
        owned=effective(owner)
        do receiver=1,size(flux)
          if(tree%topology%nodes(receiver)%level/=depth.or.origin(seam,receiver)/=0) cycle
          offset=origin(:,receiver)-origin(:,owner)
          offset(seam)=offset(seam)+domain(seam,depth+1)
          if(.not.seam_overlap(offset,tree%topology%nodes(owner)%shape, &
            tree%topology%nodes(receiver)%shape,seam)) cycle
          call reconcile_reactive_eb_peer_flux_3d(species,tree%nodes(owner)%geometry,tree%nodes(receiver)%geometry, &
            int(offset),owned,effective(receiver),ok=valid,periodic_length=period, &
        component_count=nv)
          if(.not.valid) return
        end do
      end do
    end do
    do depth=last_level,first_level,-1
    do node=1,size(flux)
      if (tree%topology%nodes(node)%level/=depth) cycle
      d=tree%topology%nodes(node)%shape
      allocate(refined(d(1),d(2),d(3)));refined=.false.
      do child=node+1,size(flux)
        if (tree%topology%nodes(child)%parent/=node) cycle
        call bounds(tree%topology,child,lo,hi)
        refined(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.true.
      end do
      if(present(reflux_increment).or.present(rereflux_delta)) then
        call compute_reactive_eb_reflux_increment_3d(tree%nodes(node)%geometry,refined,flux(node),effective(node), &
          candidate_increment%fields(node)%values,valid,all_faces=.true.)
      else
      call reflux_reactive_eb_state_3d(species,tree%nodes(node)%geometry,refined,flux(node),effective(node), &
        candidate%nodes(node)%state,candidate%nodes(node)%temperature,valid,all_faces=.true.,frame=frame)
      endif
      if (.not.valid) return
      deallocate(refined)
      if (depth==first_level) cycle
      parent=tree%topology%nodes(node)%parent;ratio=tree%topology%nodes(node)%patch%refinement_ratio
      call bounds(tree%topology,node,lo,hi)
      call restrict_eb_integrated_flux_3d(tree%nodes(node)%geometry,ratio,effective(node),restricted,valid)
      if (.not.valid) return
      effective(parent)%x(:,lo(1)-1:hi(1),lo(2):hi(2),lo(3):hi(3))=restricted%x
      effective(parent)%y(:,lo(1):hi(1),lo(2)-1:hi(2),lo(3):hi(3))=restricted%y
      effective(parent)%z(:,lo(1):hi(1),lo(2):hi(2),lo(3)-1:hi(3))=restricted%z
      call coarsen_eb_geometry_3d(tree%nodes(node)%geometry,ratio,section,valid)
      if (.not.valid) return
      section_origin=origin(:,node)/int(ratio,int64)
      do peer=1,size(flux)
        if (peer==parent.or.tree%topology%nodes(peer)%level/=depth-1) cycle
        offset=origin(:,peer)-section_origin
        if (any(offset < -int(tree%topology%nodes(peer)%shape,int64)).or.any(offset>int(hi-lo+1,int64))) cycle
        call reconcile_reactive_eb_peer_flux_3d(species,section,tree%nodes(peer)%geometry,int(offset), &
          restricted,effective(peer),ok=valid, &
        component_count=nv)
        if (.not.valid) return
      end do
      ! A restricted fine face also replaces the opposite periodic parent
      ! face, including a root that is both parent and opposite receiver.
      do seam=1,3
        if(faces(2*seam-1)/='periodic') cycle
        period=0;period(seam)=length(seam)
        do side=0,1
          if(side==0.and.section_origin(seam)/=0) cycle
          if(side==1.and.section_origin(seam)+hi(seam)-lo(seam)+1/=domain(seam,depth)) cycle
          do peer=1,size(flux)
            if(tree%topology%nodes(peer)%level/=depth-1) cycle
            if(side==0.and.origin(seam,peer)+tree%topology%nodes(peer)%shape(seam)/=domain(seam,depth)) cycle
            if(side==1.and.origin(seam,peer)/=0) cycle
            offset=origin(:,peer)-section_origin
            offset(seam)=offset(seam)+int(2*side-1,int64)*domain(seam,depth)
            if(.not.seam_overlap(offset,hi-lo+1,tree%topology%nodes(peer)%shape,seam)) cycle
            call reconcile_reactive_eb_peer_flux_3d(species,section,tree%nodes(peer)%geometry,int(offset), &
              restricted,effective(peer),ok=valid,periodic_length=period, &
        component_count=nv)
            if(.not.valid) return
          end do
        end do
      end do
    end do
    end do
    if(present(rereflux_delta)) then
      do node=1,size(tree%nodes)
        d=tree%topology%nodes(node)%shape
        do k=1,d(3)
        do j=1,d(2)
        do i=1,d(1)
          call add_eb_tree_delta_cell_3d(candidate_increment%fields(node)%values(:,i,j,k), &
            routed_delta%fields(node)%values(:,i,j,k),valid)
          if(.not.valid) return
        enddo
        enddo
        enddo
      enddo
      if(.not.present(reflux_increment)) then
        call apply_eb_tree_increment_3d(species,candidate,candidate_increment,valid,frame)
        if(.not.valid) return
      endif
    endif
    if(present(reflux_increment)) then
      call move_alloc(candidate_increment%fields,reflux_increment%fields)
      flux=effective;ok=.true.
      return
    endif
    call average_down_reactive_eb_tree_3d(species,candidate,valid,frame)
    if (.not.valid) return
    tree=candidate;flux=effective;ok=.true.
  contains
    pure logical function seam_overlap(displacement,owner_shape,receiver_shape,normal) result(overlap)
      integer(int64), intent(in) :: displacement(3)
      integer, intent(in) :: owner_shape(3),receiver_shape(3),normal
      integer :: tangent
      overlap=.true.
      do tangent=1,3
        if(tangent==normal) cycle
        if(displacement(tangent)>=owner_shape(tangent).or. &
          displacement(tangent)+receiver_shape(tangent)<=0) overlap=.false.
      end do
    end function seam_overlap
  end subroutine synchronize_reactive_eb_tree_flux_core_3d

  ! Synchronize one independent conserved density using the same EB face
  ! ownership, restriction and re-reflux ordering as the gas register.
  subroutine synchronize_reactive_eb_tree_scalar_flux_3d(species,tree,field,flux,ok, &
      minimum_level,maximum_level,boundary_faces,rereflux_delta,rereflux_first_level)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(inout) :: field
    type(eb_integrated_flux_3d), intent(inout) :: flux(:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: minimum_level,maximum_level
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(amr_tree_state_3d), optional, intent(in) :: rereflux_delta
    integer, optional, intent(in) :: rereflux_first_level
    type(amr_tree_state_3d) :: candidate,routed
    type(eb_integrated_flux_3d), allocatable :: effective(:)
    type(eb_integrated_flux_3d) :: restricted,owned
    type(eb_geometry_3d) :: section
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer(int64) :: offset(3),section_origin(3)
    real(dp), allocatable :: delta(:,:,:,:)
    real(dp) :: period(3),length(3)
    logical, allocatable :: refined(:,:,:)
    character(len=16) :: faces(6)
    integer :: p,q,owner,receiver,axis,depth,parent,child,lo(3),hi(3),d(3),ratio, &
      first_level,last_level,register_first_level,seam,side,i,j,k
    logical :: valid
    ok=.false.
    if(.not.valid_tree(tree,reactive_nvar(size(species)))) return
    if(.not.field%is_valid(tree%topology)) return
    if(size(field%fields(1)%values,1)/=1.or.size(flux)/=size(tree%nodes)) return
    faces='outflow'
    if(present(boundary_faces)) faces=boundary_faces
    if(any(faces/='outflow'.and.faces/='inflow'.and.faces/='periodic'.and. &
      faces/='slip_wall'.and.faces/='no_slip_wall'.and.faces/='isothermal_wall')) return
    if(.not.periodic_eb_faces_match(tree%nodes(1)%geometry,faces)) return
    first_level=0;last_level=tree%topology%level_count()-1
    if(present(minimum_level)) first_level=minimum_level
    if(present(maximum_level)) last_level=maximum_level
    if(first_level<0.or.last_level<first_level.or.last_level>=tree%topology%level_count()) return
    register_first_level=first_level
    if(present(rereflux_first_level)) register_first_level=rereflux_first_level
    if(register_first_level<0.or.register_first_level>first_level) return
    do p=1,size(flux)
      if(.not.valid_eb_integrated_flux_3d(flux(p),tree%nodes(p)%geometry,1)) return
    enddo
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if(.not.valid) return
    associate(g=>tree%nodes(1)%geometry)
    length=[g%x_upper-g%x_lower,g%y_upper-g%y_lower,g%z_upper-g%z_lower]
    end associate
    candidate=field;effective=flux
    if(present(rereflux_delta)) then
      call route_eb_tree_rereflux_3d(tree,rereflux_delta,routed,1,register_first_level,last_level,valid)
      if(.not.valid) return
    endif
    do p=1,size(flux)
    do q=p+1,size(flux)
      if(tree%topology%nodes(p)%level<first_level.or.tree%topology%nodes(p)%level>last_level) cycle
      if(tree%topology%nodes(q)%level/=tree%topology%nodes(p)%level) cycle
      owner=p;receiver=q
      do axis=1,3
        if(origin(axis,p)==origin(axis,q)) cycle
        if(origin(axis,p)>origin(axis,q)) then
          owner=q;receiver=p
        endif
        exit
      enddo
      offset=origin(:,receiver)-origin(:,owner)
      if(any(offset < -int(tree%topology%nodes(receiver)%shape,int64)).or. &
        any(offset>int(tree%topology%nodes(owner)%shape,int64))) cycle
      call reconcile_reactive_eb_peer_flux_3d(species,tree%nodes(owner)%geometry, &
        tree%nodes(receiver)%geometry,int(offset),effective(owner),effective(receiver),ok=valid,component_count=1)
      if(.not.valid) return
    enddo
    enddo
    do seam=1,3
      if(faces(2*seam-1)/='periodic') cycle
      period=0;period(seam)=length(seam)
      do owner=1,size(flux)
        depth=tree%topology%nodes(owner)%level
        if(depth<first_level.or.depth>last_level) cycle
        if(origin(seam,owner)+tree%topology%nodes(owner)%shape(seam)/=domain(seam,depth+1)) cycle
        owned=effective(owner)
        do receiver=1,size(flux)
          if(tree%topology%nodes(receiver)%level/=depth.or.origin(seam,receiver)/=0) cycle
          offset=origin(:,receiver)-origin(:,owner)
          offset(seam)=offset(seam)+domain(seam,depth+1)
          if(.not.seam_overlap(offset,tree%topology%nodes(owner)%shape, &
            tree%topology%nodes(receiver)%shape,seam)) cycle
          call reconcile_reactive_eb_peer_flux_3d(species,tree%nodes(owner)%geometry, &
            tree%nodes(receiver)%geometry,int(offset),owned,effective(receiver),ok=valid, &
            periodic_length=period,component_count=1)
          if(.not.valid) return
        enddo
      enddo
    enddo
    do depth=last_level,first_level,-1
    do p=1,size(flux)
      if(tree%topology%nodes(p)%level/=depth) cycle
      d=tree%topology%nodes(p)%shape
      allocate(refined(d(1),d(2),d(3)),source=.false.)
      do child=p+1,size(flux)
        if(tree%topology%nodes(child)%parent/=p) cycle
        call bounds(tree%topology,child,lo,hi)
        refined(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.true.
      enddo
      call compute_reactive_eb_reflux_increment_3d(tree%nodes(p)%geometry,refined,flux(p),effective(p), &
        delta,valid,all_faces=.true.)
      if(.not.valid) return
      if(.not.all(ieee_is_finite(delta))) return
      do k=1,d(3)
      do j=1,d(2)
      do i=1,d(1)
        if(refined(i,j,k).or.tree%nodes(p)%geometry%volume_fraction(i,j,k)<=0) cycle
        call add_eb_tree_delta_cell_3d(candidate%fields(p)%values(:,i,j,k),delta(:,i,j,k),valid)
        if(.not.valid) return
      enddo
      enddo
      enddo
      deallocate(refined,delta)
      if(depth==first_level) cycle
      parent=tree%topology%nodes(p)%parent;ratio=tree%topology%nodes(p)%patch%refinement_ratio
      call bounds(tree%topology,p,lo,hi)
      call restrict_eb_integrated_flux_3d(tree%nodes(p)%geometry,ratio,effective(p),restricted,valid)
      if(.not.valid) return
      effective(parent)%x(:,lo(1)-1:hi(1),lo(2):hi(2),lo(3):hi(3))=restricted%x
      effective(parent)%y(:,lo(1):hi(1),lo(2)-1:hi(2),lo(3):hi(3))=restricted%y
      effective(parent)%z(:,lo(1):hi(1),lo(2):hi(2),lo(3)-1:hi(3))=restricted%z
      call coarsen_eb_geometry_3d(tree%nodes(p)%geometry,ratio,section,valid)
      if(.not.valid) return
      section_origin=origin(:,p)/int(ratio,int64)
      do q=1,size(flux)
        if(q==parent.or.tree%topology%nodes(q)%level/=depth-1) cycle
        offset=origin(:,q)-section_origin
        if(any(offset < -int(tree%topology%nodes(q)%shape,int64)).or.any(offset>int(hi-lo+1,int64))) cycle
        call reconcile_reactive_eb_peer_flux_3d(species,section,tree%nodes(q)%geometry,int(offset), &
          restricted,effective(q),ok=valid,component_count=1)
        if(.not.valid) return
      enddo
      do seam=1,3
        if(faces(2*seam-1)/='periodic') cycle
        period=0;period(seam)=length(seam)
        do side=0,1
          if(side==0.and.section_origin(seam)/=0) cycle
          if(side==1.and.section_origin(seam)+hi(seam)-lo(seam)+1/=domain(seam,depth)) cycle
          do q=1,size(flux)
            if(tree%topology%nodes(q)%level/=depth-1) cycle
            if(side==0.and.origin(seam,q)+tree%topology%nodes(q)%shape(seam)/=domain(seam,depth)) cycle
            if(side==1.and.origin(seam,q)/=0) cycle
            offset=origin(:,q)-section_origin
            offset(seam)=offset(seam)+int(2*side-1,int64)*domain(seam,depth)
            if(.not.seam_overlap(offset,hi-lo+1,tree%topology%nodes(q)%shape,seam)) cycle
            call reconcile_reactive_eb_peer_flux_3d(species,section,tree%nodes(q)%geometry,int(offset), &
              restricted,effective(q),ok=valid,periodic_length=period,component_count=1)
            if(.not.valid) return
          enddo
        enddo
      enddo
    enddo
    enddo
    if(present(rereflux_delta)) then
      do p=1,size(tree%nodes)
        d=tree%topology%nodes(p)%shape
        do k=1,d(3)
        do j=1,d(2)
        do i=1,d(1)
          if(tree%nodes(p)%geometry%volume_fraction(i,j,k)<=0) cycle
          call add_eb_tree_delta_cell_3d(candidate%fields(p)%values(:,i,j,k), &
            routed%fields(p)%values(:,i,j,k),valid)
          if(.not.valid) return
        enddo
        enddo
        enddo
      enddo
    endif
    call average_down_reactive_eb_tree_field_3d(tree,candidate,valid)
    if(.not.valid) return
    field=candidate;flux=effective;ok=.true.
  contains
    pure logical function seam_overlap(displacement,owner_shape,receiver_shape,normal) result(overlap)
      integer(int64), intent(in) :: displacement(3)
      integer, intent(in) :: owner_shape(3),receiver_shape(3),normal
      integer :: tangent
      overlap=.true.
      do tangent=1,3
        if(tangent==normal) cycle
        if(displacement(tangent)>=owner_shape(tangent).or. &
          displacement(tangent)+receiver_shape(tangent)<=0) overlap=.false.
      enddo
    end function seam_overlap
  end subroutine synchronize_reactive_eb_tree_scalar_flux_3d

  subroutine synchronize_reactive_eb_tree_scalar_register_flux_3d(species,tree,field,flux, &
      extended_geometry,registers,depth,ok,boundary_faces)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(inout) :: field
    type(eb_integrated_flux_3d), intent(inout) :: flux(:)
    type(eb_geometry_3d), intent(in) :: extended_geometry(:)
    type(state_redist_register_3d), intent(inout) :: registers(:)
    integer, intent(in) :: depth
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(amr_tree_state_3d) :: correction,spill,candidate
    type(eb_integrated_flux_3d), allocatable :: working(:)
    integer :: p,first_depth
    logical :: has_spill
    ok=.false.
    if(depth<0.or.depth>=tree%topology%level_count()-1) return
    call assemble_state_redist_interface_rereflux_3d(tree%topology,extended_geometry,registers,1,depth, &
      correction,ok,boundary_faces,spill,has_spill)
    if(.not.ok) return
    first_depth=depth
    if(has_spill) then
      first_depth=depth-1
      do p=1,size(tree%nodes)
        if(tree%topology%nodes(p)%level==first_depth) correction%fields(p)%values=spill%fields(p)%values
      enddo
    endif
    candidate=field;working=flux
    call synchronize_reactive_eb_tree_scalar_flux_3d(species,tree,candidate,working,ok,depth,depth+1, &
      boundary_faces,correction,rereflux_first_level=first_depth)
    if(.not.ok) return
    field=candidate;flux=working
    do p=1,size(tree%nodes)
      if(tree%topology%nodes(p)%level==depth) then
        registers(p)%coarse_delta=0
      else if(tree%topology%nodes(p)%level==depth+1) then
        registers(p)%fine_delta=0
      endif
    enddo
  end subroutine synchronize_reactive_eb_tree_scalar_register_flux_3d

  ! Route the covered portion of an already assembled coarse correction to
  ! each child as a density increment. Parent-before-child ordering propagates
  ! it through arbitrarily nested covered cells without counting them twice.
  subroutine route_eb_tree_rereflux_3d(tree,source,increment,nv,first_level,last_level,ok)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: source
    type(amr_tree_state_3d), intent(inout) :: increment
    integer, intent(in) :: nv,first_level,last_level
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: candidate
    logical, allocatable :: inside(:,:,:)
    integer :: node,child,d(3),lo(3),hi(3),ratio
    logical :: valid
    ok=.false.
    if(.not.allocated(source%fields)) return
    if(size(source%fields)/=size(tree%nodes)) return
    do node=1,size(tree%nodes)
      d=tree%topology%nodes(node)%shape
      if(.not.allocated(source%fields(node)%values)) return
      if(any(shape(source%fields(node)%values)/=[nv,d])) return
      if(any(.not.ieee_is_finite(source%fields(node)%values))) return
      if(tree%topology%nodes(node)%level<first_level.or.tree%topology%nodes(node)%level>last_level) then
        if(any(source%fields(node)%values/=0)) return
      endif
    enddo
    candidate=source
    do node=1,size(tree%nodes)
      do child=node+1,size(tree%nodes)
        if(tree%topology%nodes(child)%parent/=node) cycle
        call bounds(tree%topology,child,lo,hi)
        ratio=tree%topology%nodes(child)%patch%refinement_ratio
        allocate(inside(hi(1)-lo(1)+1,hi(2)-lo(2)+1,hi(3)-lo(3)+1),source=.true.)
        call add_eb_rereflux_to_fine_3d( &
          candidate%fields(node)%values(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)),inside,[ratio,ratio,ratio], &
          [0,0,0],[0,0,0],candidate%fields(child)%values,valid)
        if(.not.valid) return
        candidate%fields(node)%values(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=0
        deallocate(inside)
      enddo
    enddo
    call move_alloc(candidate%fields,increment%fields);ok=.true.
  end subroutine route_eb_tree_rereflux_3d

  ! Called with a private tree; a rejected descendant never publishes an
  ! earlier accepted ancestor or its temperatures/effective fluxes.
  subroutine apply_eb_tree_increment_3d(species,tree,increment,ok,frame)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_reactive_tree_3d), intent(inout) :: tree
    type(amr_tree_state_3d), intent(in) :: increment
    logical, intent(out) :: ok
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp) :: primitive(reactive_nprim(size(species))),speed,temperature
    integer :: node,d(3),i,j,k
    logical :: valid
    ok=.false.
    do node=1,size(tree%nodes)
      d=tree%topology%nodes(node)%shape
      do k=1,d(3)
      do j=1,d(2)
      do i=1,d(1)
        if(tree%nodes(node)%geometry%volume_fraction(i,j,k)<=0) cycle
        if(all(increment%fields(node)%values(:,i,j,k)==0)) cycle
        call add_eb_tree_delta_cell_3d(tree%nodes(node)%state(:,i,j,k),increment%fields(node)%values(:,i,j,k),valid)
        if(.not.valid.or.any(tree%nodes(node)%state(6:,i,j,k)<0)) return
        call recover_reactive_eb_cell_3d(species,tree%nodes(node)%geometry,[i,j,k], &
          tree%nodes(node)%state(:,i,j,k),tree%nodes(node)%temperature(i,j,k),primitive,temperature,speed,valid,frame)
        if(.not.valid) return
        tree%nodes(node)%temperature(i,j,k)=temperature
      enddo
      enddo
      enddo
    enddo
    ok=.true.
  end subroutine apply_eb_tree_increment_3d

  subroutine add_eb_tree_delta_cell_3d(value,delta,ok)
    real(dp), intent(inout) :: value(:)
    real(dp), intent(in) :: delta(:)
    logical, intent(out) :: ok
    integer :: c
    ok=.false.
    if(size(value)/=size(delta)) return
    if(any(.not.ieee_is_finite(value)).or.any(.not.ieee_is_finite(delta))) return
    do c=1,size(value)
      if(delta(c)>0) then
        if(value(c)>huge(1._dp)-delta(c)) return
      else if(delta(c)<0) then
        if(value(c)<-huge(1._dp)-delta(c)) return
      endif
    enddo
    value=value+delta;ok=.true.
  end subroutine add_eb_tree_delta_cell_3d

  ! Six-face Mach < 0 retains outflow; active entries require inlet temperature.
  ! Pressure zero extrapolates the interior value; positive entries prescribe it.
  ! Requesting volume fractions fills interior/edges/corners too. Ancestor
  ! fractions use PCM donor values; physical corners use the first active inlet
  ! axis (x, then y, then z). No geometric fine-cell reconstruction is implied.
  ! Wider halos require explicit target-resolution volume fractions, including
  ! edges/corners. Actual output bounds correspond to 1-width:N+width; array
  ! assignment preserves caller bounds. Width one retains the legacy face route.
  ! Signed SDC sources use PCM ancestors and scalar extension: mirrored donor
  ! positions at walls, clamped elsewhere. No EOS, momentum sign flip or restriction.
  subroutine fill_sdc_eb_tree_source_halo_3d(tree,rates,target,halo_width,ghost,ok, &
      boundary_faces,target_volume_fraction,ancestor_start,ancestor_end,ancestor_fraction)
    type(eb_reactive_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: rates
    integer, intent(in) :: target,halo_width
    real(dp), intent(inout) :: ghost(:,:,:,:)
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: target_volume_fraction(:,:,:)
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    real(dp), allocatable :: candidate(:,:,:,:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer(int64) :: global(3),mapped(3),local(3)
    integer :: boundary_face
    integer :: n(3),index(3),source(3),nv,width,level,depth,node,donor,i,j,k,axis
    real(dp) :: fraction,volume
    character(len=32) :: faces(6)
    logical :: valid
    ok=.false.;width=halo_width;nv=size(ghost,1)
    if(nv<1.or.width<1.or.width>7) return
    if(.not.tree%topology%is_valid().or..not.allocated(tree%nodes)) return
    if(size(tree%nodes)/=size(tree%topology%nodes)) return
    if(target<1.or.target>size(tree%nodes)) return
    if(.not.allocated(tree%nodes(1)%state)) return
    if(.not.valid_tree(tree,size(tree%nodes(1)%state,1))) return
    if(.not.field_shape(rates)) return
    faces='outflow'
    if(present(boundary_faces)) faces=boundary_faces
    if(.not.valid_reactive_boundaries_3d(faces)) return
    if(.not.periodic_eb_faces_match(tree%nodes(1)%geometry,faces)) return
    n=tree%topology%nodes(target)%shape;level=tree%topology%nodes(target)%level+1
    if(any(n>huge(1)-2*width)) return
    if(any(shape(ghost)/=[nv,n+2*width])) return
    if(present(target_volume_fraction)) then
      if(any(shape(target_volume_fraction)/=n+2*width)) return
      if(.not.all(ieee_is_finite(target_volume_fraction))) return
      if(any(target_volume_fraction<0).or.any(target_volume_fraction>1)) return
      if(any(target_volume_fraction(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))/= &
        tree%nodes(target)%geometry%volume_fraction)) return
    end if
    if(present(ancestor_start).or.present(ancestor_end).or.present(ancestor_fraction)) then
      if(.not.(present(ancestor_start).and.present(ancestor_end).and.present(ancestor_fraction))) return
      if(.not.field_shape(ancestor_start).or..not.field_shape(ancestor_end)) return
      if(size(ancestor_fraction)/=level-1) return
      if(.not.all(ieee_is_finite(ancestor_fraction))) return
      if(any(ancestor_fraction<0).or.any(ancestor_fraction>1)) return
    end if
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if(.not.valid) return
    allocate(candidate(nv,n(1)+2*width,n(2)+2*width,n(3)+2*width),source=0._dp)
    do k=1,n(3)+2*width
    do j=1,n(2)+2*width
    do i=1,n(1)+2*width
      index=[i,j,k]-width
      global=origin(:,target)+int(index-1,int64)
      do axis=1,3
        if(faces(2*axis-1)=='periodic') global(axis)=modulo(global(axis),domain(axis,level))
      end do
      ! Match state/geometry donor locations at physical walls. Source
      ! components retain their existing scalar extension convention.
      do axis=1,3
        do while(global(axis)<0.or.global(axis)>=domain(axis,level))
          boundary_face=2*axis
          if(global(axis)<0) boundary_face=boundary_face-1
          if(faces(boundary_face)/='slip_wall'.and.faces(boundary_face)/='no_slip_wall'.and. &
            faces(boundary_face)/='isothermal_wall') exit
          if(global(axis)<0) then
            global(axis)=-global(axis)-1
          else
            global(axis)=2*domain(axis,level)-global(axis)-1
          end if
        end do
      end do
      global=max(0_int64,min(domain(:,level)-1,global))
      donor=0
      if(all(index>=1).and.all(index<=n)) then
        donor=target;source=index
      else
        do depth=level,1,-1
          mapped=global/(domain(:,level)/domain(:,depth))
          do node=1,size(tree%nodes)
            if(tree%topology%nodes(node)%level+1/=depth) cycle
            local=mapped-origin(:,node)
            if(any(local<0).or.any(local>=int(tree%topology%nodes(node)%shape,int64))) cycle
            donor=node;source=int(local)+1;exit
          end do
          if(donor/=0) exit
        end do
      end if
      if(donor==0) return
      volume=tree%nodes(donor)%geometry%volume_fraction(source(1),source(2),source(3))
      if(present(target_volume_fraction)) then
        if(target_volume_fraction(i,j,k)<=0) cycle
        if(volume<=0) return
      else
        if(volume<=0) cycle
      end if
      depth=tree%topology%nodes(donor)%level+1
      if(present(ancestor_fraction).and.depth<level) then
        fraction=ancestor_fraction(depth)
        if(fraction==0) then
          candidate(:,i,j,k)=ancestor_start%fields(donor)%values(:,source(1),source(2),source(3))
        else if(fraction==1) then
          candidate(:,i,j,k)=ancestor_end%fields(donor)%values(:,source(1),source(2),source(3))
        else
          candidate(:,i,j,k)=(1-fraction)*ancestor_start%fields(donor)%values(:,source(1),source(2),source(3))+ &
            fraction*ancestor_end%fields(donor)%values(:,source(1),source(2),source(3))
        end if
      else
        candidate(:,i,j,k)=rates%fields(donor)%values(:,source(1),source(2),source(3))
      end if
      if(.not.all(ieee_is_finite(candidate(:,i,j,k)))) return
    end do
    end do
    end do
    ghost=candidate;ok=.true.
  contains
    logical function field_shape(field) result(good)
      type(amr_tree_state_3d), intent(in) :: field
      integer :: p
      good=.false.
      if(.not.allocated(field%fields)) return
      if(size(field%fields)/=size(tree%nodes)) return
      do p=1,size(field%fields)
        if(.not.allocated(field%fields(p)%values)) return
        if(any(shape(field%fields(p)%values)/=[nv,tree%topology%nodes(p)%shape])) return
      end do
      good=.true.
    end function field_shape
  end subroutine fill_sdc_eb_tree_source_halo_3d

  subroutine fill_reactive_eb_tree_halo_3d(species,tree,target,ghost,ghost_temperature,ok, &
      ancestor_start,ancestor_end,ancestor_fraction,inflow_mach,inflow_temperature,inflow_pressure, &
      ghost_volume_fraction,halo_width,target_volume_fraction,boundary_faces,inflow_state,restrict_snapshot,frame, &
      internal_energy,ghost_internal_energy,ancestor_internal_energy_start,ancestor_internal_energy_end,coarse_geometry, &
      separate_momentum,consistent_velocity,paired_pressure,pressure_level,connected_directions,preserve_temperature_seed)
    logical, optional, intent(in) :: preserve_temperature_seed
    logical :: keep_seed
    integer, optional, intent(in) :: pressure_level
    type(eb_geometry_3d), optional, intent(in) :: coarse_geometry(:)
    logical, optional, intent(in) :: separate_momentum,consistent_velocity,paired_pressure,connected_directions
    type(nasa7_species), intent(in) :: species(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(eb_reactive_tree_3d), intent(in) :: tree
    integer, intent(in) :: target
    integer, optional, intent(in) :: halo_width
    real(dp), optional, intent(in) :: target_volume_fraction(:,:,:)
    real(dp), intent(inout) :: ghost(:,0:,0:,0:),ghost_temperature(0:,0:,0:)
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    real(dp), optional, intent(inout) :: ghost_volume_fraction(0:,0:,0:)
    real(dp), allocatable :: volumes(:,:,:)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    character(len=16) :: faces(6),linear_faces(6)
    integer :: inlet_face,momentum_sign(3)
    logical, allocatable :: filled(:,:,:)
    logical, intent(out) :: ok
    type(eb_reactive_tree_3d), optional, intent(in) :: ancestor_start,ancestor_end
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy,ancestor_internal_energy_start, &
      ancestor_internal_energy_end
    real(dp), optional, intent(inout) :: ghost_internal_energy(0:,0:,0:)
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    type(eb_reactive_tree_3d) :: snapshot
    type(amr_tree_state_3d) :: held_snapshot,pressure_field
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: node_shapes(:,:)
    integer :: shape_node
    integer, allocatable :: ratios(:)
    integer(int64) :: global(3),mapped(3),local(3),offset(3)
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:),primitive(:),held(:,:,:),cell_held
    real(dp) :: area,guess,recovered,speed,fraction,value(reactive_nvar(size(species)))
    integer :: n(3),nv,level,depth,node,donor,index(3),inside(3),source(3),first(3),last(3),axis,side,i,j,k,width
    real(dp) :: mach(6),pin(6),tin(6),inlet_value(reactive_nvar(size(species))),inlet_t
    logical, optional, intent(in) :: restrict_snapshot
    logical :: valid,physical,restrict_state,linear,pressure_policy,restrict_pressure
    integer :: stencil(4,6),stencil_signs(3,6),sample_node,sample_index(3),q
    real(dp) :: packed(reactive_nvar(size(species))+2),neighbors(reactive_nvar(size(species))+2,6),shift(3)
    real(dp) :: target_pressure(1),pressure_neighbors(1,6),energy_delta
    integer :: pp,xx,yy,zz,mm(3)
    keep_seed=.false.
    if(present(preserve_temperature_seed)) keep_seed=preserve_temperature_seed
    ok=.false.;nv=reactive_nvar(size(species));width=1
    restrict_pressure=.not.present(ancestor_fraction)
    if(present(pressure_level)) restrict_pressure=restrict_pressure.and.pressure_level<0
    pressure_policy=.false.
    if(present(paired_pressure)) pressure_policy=paired_pressure
    if(pressure_policy) then
      if(.not.present(consistent_velocity)) return
      if(.not.consistent_velocity) return
      if(present(frame)) return
    endif
    if(present(consistent_velocity)) then
      if(consistent_velocity) then
        if(.not.present(separate_momentum)) return
        if(.not.separate_momentum) return
      endif
    endif
    if(present(connected_directions)) then
      if(connected_directions) then
        if(.not.present(separate_momentum)) return
        if(.not.separate_momentum) return
      endif
    endif
    if(present(separate_momentum)) then
      if(separate_momentum.and..not.present(coarse_geometry)) return
    endif
    if(present(frame)) then
      if(.not.frame%is_valid()) return
    end if
    if(present(halo_width)) width=halo_width
    if(width<1.or.width>huge(1)/4) return
    if(width>1.and.(.not.present(ghost_volume_fraction).or..not.present(target_volume_fraction))) return
    if(present(internal_energy).neqv.present(ghost_internal_energy)) return
    if(present(internal_energy)) then
      if(present(ancestor_internal_energy_start).neqv.present(ancestor_internal_energy_end)) return
      if(present(ancestor_fraction).neqv.present(ancestor_internal_energy_start)) return
    else
      if(present(ancestor_internal_energy_start).or.present(ancestor_internal_energy_end)) return
    endif
    if(present(coarse_geometry)) then
      if(.not.present(internal_energy).or.width<=1) return
    endif
    faces='outflow'
    if(present(boundary_faces)) then
      if(any(boundary_faces/='inflow'.and.boundary_faces/='outflow'.and.boundary_faces/='periodic'.and. &
        boundary_faces/='slip_wall'.and.boundary_faces/='no_slip_wall'.and.boundary_faces/='isothermal_wall')) return
      if(present(inflow_mach).or.present(inflow_pressure)) return
      faces=boundary_faces
      if(any(faces=='inflow')) then
        if(.not.present(inflow_state).or..not.present(inflow_temperature)) return
      end if
    else
      if(present(inflow_state)) return
    end if
    if(present(inflow_state)) then
      if(any(shape(inflow_state)/=[nv,6])) return
      if(.not.all(ieee_is_finite(inflow_state))) return
    end if
    mach=-1;pin=0;tin=0
    if (present(inflow_temperature).or.present(inflow_pressure)) then
      if (.not.present(inflow_mach).and..not.present(boundary_faces)) return
    end if
    if (present(inflow_mach)) mach=inflow_mach
    if (present(inflow_temperature)) tin=inflow_temperature
    if (present(inflow_pressure)) pin=inflow_pressure
    if (.not.all(ieee_is_finite(mach)).or.any(mach>100)) return
    if (.not.all(ieee_is_finite(pin)).or.any(pin<0).or.any(pin>1e15_dp)) return
    if (.not.all(ieee_is_finite(tin))) return
    if (any(pin>0.and.mach<0).or.any(tin<=0.and.mach>=0)) return
    if (.not.valid_tree(tree,nv)) return
    if(present(coarse_geometry)) then
      if(size(coarse_geometry)/=size(tree%nodes)) return
    endif
    if(.not.periodic_eb_faces_match(tree%nodes(1)%geometry,faces)) return
    if (target<1.or.target>size(tree%nodes)) return
    n=tree%topology%nodes(target)%shape;level=tree%topology%nodes(target)%level+1
    if(present(internal_energy)) then
      if(.not.valid_held(internal_energy)) return
      if(any(shape(ghost_internal_energy)/=n+2*width)) return
      if(present(ancestor_internal_energy_start)) then
        if(.not.valid_held(ancestor_internal_energy_start).or. &
          .not.valid_held(ancestor_internal_energy_end)) return
      endif
    endif
    if(any(n>huge(1)-2*width)) return
    if (any(shape(ghost)/=[nv,n+2*width]).or.any(shape(ghost_temperature)/=n+2*width)) return
    if(present(ghost_volume_fraction)) then
      if(any(shape(ghost_volume_fraction)/=n+2*width)) return
    end if
    if(present(target_volume_fraction)) then
      if(.not.present(ghost_volume_fraction)) return
      if(any(shape(target_volume_fraction)/=n+2*width)) return
      if(.not.all(ieee_is_finite(target_volume_fraction))) return
      if(any(target_volume_fraction<0).or.any(target_volume_fraction>1)) return
      if(any(target_volume_fraction(width+1:width+n(1),width+1:width+n(2),width+1:width+n(3))/= &
        tree%nodes(target)%geometry%volume_fraction)) return
    end if
    linear_faces=faces
    where(mach>=0) linear_faces='inflow'
    allocate(node_shapes(3,size(tree%topology%nodes)))
    do shape_node=1,size(tree%topology%nodes)
      node_shapes(:,shape_node)=tree%topology%nodes(shape_node)%shape
    enddo
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if (.not.valid) return
    if (present(ancestor_start).or.present(ancestor_end).or.present(ancestor_fraction)) then
      if (.not.(present(ancestor_start).and.present(ancestor_end).and.present(ancestor_fraction))) return
      if (size(ancestor_fraction)/=level-1) return
      if (.not.all(ieee_is_finite(ancestor_fraction))) return
      if (any(ancestor_fraction<0._dp).or.any(ancestor_fraction>1._dp)) return
      if (.not.matching_snapshot(ancestor_start).or..not.matching_snapshot(ancestor_end)) return
    end if
    snapshot=tree
    restrict_state=.true.
    if(present(restrict_snapshot)) restrict_state=restrict_snapshot
    if(present(internal_energy)) then
      held_snapshot=internal_energy
      if(restrict_state) then
        call average_down_reactive_eb_tree_3d(species,snapshot,valid,frame,held_snapshot)
        if(.not.valid) return
      endif
    else if(restrict_state) then
      call average_down_reactive_eb_tree_3d(species,snapshot,valid,frame)
      if(.not.valid) return
    endif
    if(pressure_policy) then
      ! Evaluate EOS before restriction: p(mean U) differs from mean p(U).
      pressure_field=internal_energy
      allocate(primitive(reactive_nprim(size(species))))
      do pp=1,size(tree%nodes)
        mm=tree%topology%nodes(pp)%shape
        pressure_field%fields(pp)%values=0
        do zz=1,mm(3)
        do yy=1,mm(2)
        do xx=1,mm(1)
          if(tree%nodes(pp)%geometry%volume_fraction(xx,yy,zz)<=0) cycle
          call reactive_conserved_to_primitive(species,tree%nodes(pp)%state(:,xx,yy,zz), &
            tree%nodes(pp)%temperature(xx,yy,zz),primitive,recovered,speed,valid, &
            internal_energy_density=internal_energy%fields(pp)%values(1,xx,yy,zz))
          if(.not.valid) return
          pressure_field%fields(pp)%values(1,xx,yy,zz)=primitive(5)
        enddo
        enddo
        enddo
      enddo
      deallocate(primitive)
      if(restrict_pressure) then
        call average_down_reactive_eb_tree_field_3d(tree,pressure_field,valid)
        if(.not.valid) return
      endif
    endif
    allocate(candidate(nv,1-width:n(1)+width,1-width:n(2)+width,1-width:n(3)+width), &
      temps(1-width:n(1)+width,1-width:n(2)+width,1-width:n(3)+width), &
      primitive(reactive_nprim(size(species))))
    candidate=0._dp;temps=0._dp
    if(present(internal_energy)) then
      allocate(held(1-width:n(1)+width,1-width:n(2)+width,1-width:n(3)+width),cell_held)
      held=0._dp
    endif
    allocate(filled(1-width:n(1)+width,1-width:n(2)+width,1-width:n(3)+width));filled=.false.
    if(width==1) then
    do axis=1,3
    do side=0,1
      first=1;last=n;first(axis)=side*(n(axis)+1);last(axis)=first(axis)
      do k=first(3),last(3)
      do j=first(2),last(2)
      do i=first(1),last(1)
        index=[i,j,k];inside=max(1,min(n,index))
        select case(axis)
        case(1);area=tree%nodes(target)%geometry%x_face_fraction(side*n(1),j,k)
        case(2);area=tree%nodes(target)%geometry%y_face_fraction(i,side*n(2),k)
        case(3);area=tree%nodes(target)%geometry%z_face_fraction(i,j,side*n(3))
        end select
        if (area<=0._dp) cycle
        global=origin(:,target)+int(index-1,int64)
        donor=0
        if(faces(2*axis-1)=='periodic') global(axis)=modulo(global(axis),domain(axis,level))
        physical=any(global<0_int64).or.any(global>=domain(:,level))
        if (physical) then
          ! Physical boundary starts from the target's adjacent interior cell.
          donor=target;source=inside
        else
          do depth=level,1,-1
            mapped=global/(domain(:,level)/domain(:,depth))
            do node=1,size(tree%nodes)
              if (tree%topology%nodes(node)%level+1/=depth) cycle
              local=mapped-origin(:,node)
              if (any(local<0_int64).or.any(local>=int(tree%topology%nodes(node)%shape,int64))) cycle
              donor=node;source=int(local)+1;exit
            end do
            if (donor/=0) exit
          end do
        end if
        if (donor==0) return
        if (snapshot%nodes(donor)%geometry%volume_fraction(source(1),source(2),source(3))<=0._dp) return
        guess=snapshot%nodes(donor)%temperature(source(1),source(2),source(3))
        value=snapshot%nodes(donor)%state(:,source(1),source(2),source(3))
        if(present(internal_energy)) cell_held=held_snapshot%fields(donor)%values(1,source(1),source(2),source(3))
        depth=tree%topology%nodes(donor)%level+1
        if (present(ancestor_fraction)) then
          if (depth<level) then
            fraction=ancestor_fraction(depth)
            ! Endpoint states are provisional level updates. Do not average
            ! them down from descendants that may be at a different time.
            if (fraction==0._dp) then
              value=ancestor_start%nodes(donor)%state(:,source(1),source(2),source(3))
              guess=ancestor_start%nodes(donor)%temperature(source(1),source(2),source(3))
            else if (fraction==1._dp) then
              value=ancestor_end%nodes(donor)%state(:,source(1),source(2),source(3))
              guess=ancestor_end%nodes(donor)%temperature(source(1),source(2),source(3))
            else
              value=(1._dp-fraction)*ancestor_start%nodes(donor)%state(:,source(1),source(2),source(3))+ &
                fraction*ancestor_end%nodes(donor)%state(:,source(1),source(2),source(3))
              guess=(1._dp-fraction)*ancestor_start%nodes(donor)%temperature(source(1),source(2),source(3))+ &
                fraction*ancestor_end%nodes(donor)%temperature(source(1),source(2),source(3))
            end if
            if(present(internal_energy)) cell_held= &
              (1._dp-fraction)*ancestor_internal_energy_start%fields(donor)%values(1,source(1),source(2),source(3))+ &
              fraction*ancestor_internal_energy_end%fields(donor)%values(1,source(1),source(2),source(3))
          end if
        end if
        if (.not.ieee_is_finite(guess)) return
        if (guess<=0._dp) return
        call recover_reactive_eb_cell_3d(species,tree%nodes(target)%geometry,[i,j,k],value, &
          guess,primitive,recovered,speed,valid,frame,cell_held)
        if (.not.valid) return
        if (physical.and.mach(2*axis-1+side)>=0) then
          call boundary_energy(.false.,valid)
          if(.not.valid) return
          if (pin(2*axis-1+side)>0) then
            call mach_inlet_ghost_3d(species,value,recovered,2*axis-1+side,tin(2*axis-1+side), &
              mach(2*axis-1+side),inlet_value,inlet_t,valid,inlet_pressure=pin(2*axis-1+side), &
              internal_energy_density=cell_held)
          else
            call mach_inlet_ghost_3d(species,value,recovered,2*axis-1+side,tin(2*axis-1+side), &
              mach(2*axis-1+side),inlet_value,inlet_t,valid,internal_energy_density=cell_held)
          end if
          if (.not.valid) return
          value=inlet_value;recovered=inlet_t
          call boundary_energy(.true.,valid)
          if(.not.valid) return
        end if
        if(physical) then
          call prescribed_inlet(2*axis-1+side,valid)
          if(.not.valid) return
          if(faces(2*axis-1+side)=='slip_wall') value(1+axis)=-value(1+axis)
          if(faces(2*axis-1+side)=='no_slip_wall'.or.faces(2*axis-1+side)=='isothermal_wall') value(2:4)=-value(2:4)
        end if
        if(keep_seed.and.depth<level.and..not.physical.and..not.pressure_policy) recovered=guess
        candidate(:,i,j,k)=value
        temps(i,j,k)=recovered
        if(present(internal_energy)) then
          if(physical) then
            if(mach(2*axis-1+side)>=0.or.faces(2*axis-1+side)=='inflow') then
              call gas_internal_energy(value,cell_held,valid)
              if(.not.valid) return
            endif
          endif
          held(i,j,k)=cell_held
        endif
        filled(i,j,k)=.true.
      end do
      end do
      end do
    end do
    end do
    end if
    if(present(ghost_volume_fraction)) then
      allocate(volumes(1-width:n(1)+width,1-width:n(2)+width,1-width:n(3)+width));volumes=0
      do k=1-width,n(3)+width
      do j=1-width,n(2)+width
      do i=1-width,n(1)+width
        index=[i,j,k];global=origin(:,target)+int(index-1,int64)
        inlet_face=0
        do axis=1,3
          side=-1
          if(global(axis)<0) side=0
          if(global(axis)>=domain(axis,level)) side=1
          if(side<0) cycle
          if((mach(2*axis-1+side)>=0.or.faces(2*axis-1+side)=='inflow').and.inlet_face==0) inlet_face=2*axis-1+side
        end do
        call map_wall_donor(global,momentum_sign)
        global=max(0_int64,min(domain(:,level)-1,global))
        donor=0
        if(all(index>=1).and.all(index<=n)) then
          donor=target;source=index
        else
          do depth=level,1,-1
            mapped=global/(domain(:,level)/domain(:,depth))
            do node=1,size(tree%nodes)
              if(tree%topology%nodes(node)%level+1/=depth) cycle
              local=mapped-origin(:,node)
              if(any(local<0).or.any(local>=int(tree%topology%nodes(node)%shape,int64))) cycle
              donor=node;source=int(local)+1;exit
            end do
            if(donor/=0) exit
          end do
        end if
        if(donor==0) return
        volumes(i,j,k)=snapshot%nodes(donor)%geometry%volume_fraction(source(1),source(2),source(3))
        if(present(target_volume_fraction)) then
          volumes(i,j,k)=target_volume_fraction(i+width,j+width,k+width)
        else if(allocated(snapshot%nodes(target)%geometry%halo_volume_fraction)) then
          local=global-origin(:,target)+1
          do axis=1,3
            if(faces(2*axis-1)=='periodic'.or.faces(2*axis-1)=='slip_wall'.or. &
              faces(2*axis-1)=='no_slip_wall'.or.faces(2*axis)=='slip_wall'.or. &
              faces(2*axis)=='no_slip_wall'.or.any(faces(2*axis-1:2*axis)=='isothermal_wall')) local(axis)=index(axis)
          end do
          if(snapshot%nodes(target)%geometry%halo_volume_fraction(local(1),local(2),local(3))>=0) &
            volumes(i,j,k)=snapshot%nodes(target)%geometry%halo_volume_fraction(local(1),local(2),local(3))
        end if
        if(volumes(i,j,k)<=0) cycle
        if(width>1.and.snapshot%nodes(donor)%geometry%volume_fraction(source(1),source(2),source(3))<=0) return
        if(filled(i,j,k)) cycle
        value=snapshot%nodes(donor)%state(:,source(1),source(2),source(3))
        guess=snapshot%nodes(donor)%temperature(source(1),source(2),source(3))
        if(present(internal_energy)) held(i,j,k)=held_snapshot%fields(donor)%values(1,source(1),source(2),source(3))
        depth=tree%topology%nodes(donor)%level+1
        if(present(ancestor_fraction)) then
          if(depth<level) then
            fraction=ancestor_fraction(depth)
            value=(1-fraction)*ancestor_start%nodes(donor)%state(:,source(1),source(2),source(3))+ &
              fraction*ancestor_end%nodes(donor)%state(:,source(1),source(2),source(3))
            guess=(1-fraction)*ancestor_start%nodes(donor)%temperature(source(1),source(2),source(3))+ &
              fraction*ancestor_end%nodes(donor)%temperature(source(1),source(2),source(3))
            if(present(internal_energy)) held(i,j,k)= &
              (1-fraction)*ancestor_internal_energy_start%fields(donor)%values(1,source(1),source(2),source(3))+ &
              fraction*ancestor_internal_energy_end%fields(donor)%values(1,source(1),source(2),source(3))
          end if
        end if
        if(pressure_policy) then
          call sample_pressure(donor,source,target_pressure(1),valid)
          if(.not.valid) return
        endif
        if(present(coarse_geometry).and.depth<level) then
          call plan_eb_linear_halo_3d(tree%topology%nodes%level,node_shapes,origin,domain,donor,source,level,global, &
            coarse_geometry(donor),width,linear_faces,stencil,stencil_signs,shift,linear,valid, &
      connected_directions=connected_directions)
          if(.not.valid) return
          if(linear) then
            packed=[value,guess,held(i,j,k)]
            do q=1,6
              sample_node=stencil(1,q);sample_index=stencil(2:4,q)
              call linear_sample(sample_node,sample_index,neighbors(:,q))
              neighbors(2:4,q)=neighbors(2:4,q)*real(stencil_signs(:,q),dp)
              if(pressure_policy) then
                call sample_pressure(sample_node,sample_index,pressure_neighbors(1,q),valid)
                if(.not.valid) return
              endif
            enddo
            call interpolate_eb_linear_halo_3d(packed,neighbors,shift,valid,separate_momentum,consistent_velocity)
            if(.not.valid) return
            value=packed(:nv);guess=packed(nv+1);held(i,j,k)=packed(nv+2)
            if(pressure_policy) then
              call interpolate_eb_linear_halo_3d(target_pressure,pressure_neighbors,shift,valid)
              if(.not.valid) return
            endif
          endif
        endif
        if(pressure_policy) then
          if(depth<level.or.covered(tree%topology,donor,source)) then
            call reconstruct_eb_pressure_state_3d(species,value,held(i,j,k),guess,target_pressure(1),energy_delta,valid)
            if(.not.valid) return
          endif
        endif
        if(present(internal_energy)) cell_held=held(i,j,k)
        call recover_reactive_eb_cell_3d(species,tree%nodes(target)%geometry,[i,j,k],value,guess,primitive, &
          recovered,speed,valid,frame,cell_held)
        if(.not.valid) return
        if(inlet_face>0) then
          if(faces(inlet_face)=='inflow') then
            call prescribed_inlet(inlet_face,valid)
            if(.not.valid) return
          else
          call boundary_energy(.false.,valid)
          if(.not.valid) return
          if(pin(inlet_face)>0) then
            call mach_inlet_ghost_3d(species,value,recovered,inlet_face,tin(inlet_face),mach(inlet_face), &
              inlet_value,inlet_t,valid,inlet_pressure=pin(inlet_face),internal_energy_density=cell_held)
          else
            call mach_inlet_ghost_3d(species,value,recovered,inlet_face,tin(inlet_face),mach(inlet_face), &
              inlet_value,inlet_t,valid,internal_energy_density=cell_held)
          end if
          if(.not.valid) return
          value=inlet_value;recovered=inlet_t
          call boundary_energy(.true.,valid)
          if(.not.valid) return
          end if
        end if
        value(2:4)=value(2:4)*real(momentum_sign,dp)
        ! Validation must not advance the native inverse-EOS seed before
        ! the MOL primitive conversion. Inlet/pressure reconstruction owns T.
        if(keep_seed.and.depth<level.and.inlet_face==0.and..not.pressure_policy) recovered=guess
        candidate(:,i,j,k)=value;temps(i,j,k)=recovered
        if(present(internal_energy).and.inlet_face>0) then
          call gas_internal_energy(value,held(i,j,k),valid)
          if(.not.valid) return
        endif
      end do
      end do
      end do
    end if
    if(pressure_policy) then
      ! Peer overlays must use the same covered representation as the full
      ! donor pass. These copies never publish to valid leaf state.
      do pp=1,size(tree%nodes)
        mm=tree%topology%nodes(pp)%shape
        do zz=1,mm(3)
        do yy=1,mm(2)
        do xx=1,mm(1)
          if(tree%nodes(pp)%geometry%volume_fraction(xx,yy,zz)<=0) cycle
          if(.not.covered(tree%topology,pp,[xx,yy,zz])) cycle
          call reconstruct_eb_pressure_state_3d(species,snapshot%nodes(pp)%state(:,xx,yy,zz), &
            held_snapshot%fields(pp)%values(1,xx,yy,zz),snapshot%nodes(pp)%temperature(xx,yy,zz), &
            pressure_field%fields(pp)%values(1,xx,yy,zz),energy_delta,valid)
          if(.not.valid) return
        enddo
        enddo
        enddo
      enddo
    endif
    ! The peer overlay also verifies open-face aperture and centroid agreement,
    ! including peers that belong to different parent patches.
    do node=1,size(tree%nodes)
      if (node==target.or.tree%topology%nodes(node)%level+1/=level) cycle
      offset=origin(:,node)-origin(:,target)
      if (any(offset < -int(tree%topology%nodes(node)%shape,int64)).or.any(offset>int(n,int64))) cycle
      if(present(internal_energy)) then
        call overlay_reactive_eb_peer_halo_3d(species,snapshot%nodes(target)%geometry,snapshot%nodes(node)%geometry, &
          int(offset),snapshot%nodes(node)%state,snapshot%nodes(node)%temperature, &
          candidate(:,0:n(1)+1,0:n(2)+1,0:n(3)+1),temps(0:n(1)+1,0:n(2)+1,0:n(3)+1),valid,frame=frame, &
          peer_internal_energy=held_snapshot%fields(node)%values(1,:,:,:), &
          ghost_internal_energy=held(0:n(1)+1,0:n(2)+1,0:n(3)+1))
      else
        call overlay_reactive_eb_peer_halo_3d(species,snapshot%nodes(target)%geometry,snapshot%nodes(node)%geometry, &
          int(offset),snapshot%nodes(node)%state,snapshot%nodes(node)%temperature, &
          candidate(:,0:n(1)+1,0:n(2)+1,0:n(3)+1),temps(0:n(1)+1,0:n(2)+1,0:n(3)+1),valid,frame=frame)
      endif
      if (.not.valid) return
    end do
    ghost=candidate;ghost_temperature=temps
    if(present(ghost_internal_energy)) ghost_internal_energy=held
    if(present(ghost_volume_fraction)) ghost_volume_fraction=volumes
    ok=.true.
  contains
    subroutine sample_pressure(p,idx,result,good)
      integer, intent(in) :: p,idx(3)
      real(dp), intent(inout) :: result
      logical, intent(out) :: good
      integer :: x,y,z,d
      x=idx(1);y=idx(2);z=idx(3);d=tree%topology%nodes(p)%level+1
      result=pressure_field%fields(p)%values(1,x,y,z);good=.true.
      if(.not.present(ancestor_fraction)) return
      if(d>=level) return
      call interpolate_eb_pressure_time_3d(species,ancestor_start%nodes(p)%state(:,x,y,z), &
        ancestor_end%nodes(p)%state(:,x,y,z),ancestor_internal_energy_start%fields(p)%values(1,x,y,z), &
        ancestor_internal_energy_end%fields(p)%values(1,x,y,z),ancestor_start%nodes(p)%temperature(x,y,z), &
        ancestor_end%nodes(p)%temperature(x,y,z),ancestor_fraction(d),result,good)
    end subroutine sample_pressure
    subroutine linear_sample(p,idx,result)
      integer, intent(in) :: p,idx(3)
      real(dp), intent(out) :: result(nv+2)
      integer :: x,y,z
      x=idx(1);y=idx(2);z=idx(3)
      result=[snapshot%nodes(p)%state(:,x,y,z),snapshot%nodes(p)%temperature(x,y,z), &
        held_snapshot%fields(p)%values(1,x,y,z)]
      if(present(ancestor_fraction)) then
        result=(1-fraction)*[ancestor_start%nodes(p)%state(:,x,y,z),ancestor_start%nodes(p)%temperature(x,y,z), &
          ancestor_internal_energy_start%fields(p)%values(1,x,y,z)]+ &
          fraction*[ancestor_end%nodes(p)%state(:,x,y,z),ancestor_end%nodes(p)%temperature(x,y,z), &
          ancestor_internal_energy_end%fields(p)%values(1,x,y,z)]
      endif
    end subroutine linear_sample
    logical function valid_held(field) result(good)
      type(amr_tree_state_3d), intent(in) :: field
      integer :: p,m(3)
      good=.false.
      if(.not.allocated(field%fields)) return
      if(size(field%fields)/=size(tree%nodes)) return
      do p=1,size(tree%nodes)
        if(.not.allocated(field%fields(p)%values)) return
        m=tree%topology%nodes(p)%shape
        if(any(shape(field%fields(p)%values)/=[1,m])) return
        if(any(.not.ieee_is_finite(field%fields(p)%values(1,:,:,:)).and. &
          tree%nodes(p)%geometry%volume_fraction>0)) return
      enddo
      good=.true.
    end function valid_held
    subroutine gas_internal_energy(conserved,rhoe,accepted)
      real(dp), intent(in) :: conserved(:)
      real(dp), intent(out) :: rhoe
      logical, intent(out) :: accepted
      real(dp) :: cell(size(conserved),1,1,1),temp(1,1,1),volume(1,1,1),scalar(1,1,1)
      cell(:,1,1,1)=conserved;temp=recovered;volume=1;scalar=0
      call initialize_nasa7_internal_energy_patch(species,cell,temp,volume,scalar,accepted)
      rhoe=scalar(1,1,1)
    end subroutine gas_internal_energy
    ! Prescribed and Mach inlet APIs use ordinary relative total energy.
    subroutine boundary_energy(to_rotating,accepted)
      logical, intent(in) :: to_rotating
      logical, intent(out) :: accepted
      real(dp) :: position(3),energy
      accepted=.true.
      if(.not.present(frame)) return
      position=[tree%nodes(target)%geometry%x_lower,tree%nodes(target)%geometry%y_lower, &
        tree%nodes(target)%geometry%z_lower]+(real([i,j,k],dp)-.5_dp)* &
        [tree%nodes(target)%geometry%dx,tree%nodes(target)%geometry%dy,tree%nodes(target)%geometry%dz]
      energy=value(5)
      call convert_rotating_energy_3d(frame,position,value(1),value(5),to_rotating,energy,accepted)
      if(accepted) value(5)=energy
    end subroutine boundary_energy
    subroutine map_wall_donor(position,signs)
      integer(int64), intent(inout) :: position(3)
      integer, intent(out) :: signs(3)
      integer :: direction,face
      signs=1
      do direction=1,3
        if(faces(2*direction-1)=='periodic') then
          position(direction)=modulo(position(direction),domain(direction,level))
          cycle
        end if
        do while(position(direction)<0.or.position(direction)>=domain(direction,level))
          face=2*direction
          if(position(direction)<0) face=face-1
          if(faces(face)/='slip_wall'.and.faces(face)/='no_slip_wall'.and.faces(face)/='isothermal_wall') exit
          if(faces(face)=='no_slip_wall'.or.faces(face)=='isothermal_wall') then
            signs=-signs
          else
            signs(direction)=-signs(direction)
          end if
          if(position(direction)<0) then
            position(direction)=-position(direction)-1
          else
            position(direction)=2*domain(direction,level)-position(direction)-1
          end if
        end do
      end do
    end subroutine
    subroutine prescribed_inlet(face,accepted)
      integer, intent(in) :: face
      logical, intent(out) :: accepted
      accepted=.true.
      if(faces(face)/='inflow') return
      if(present(internal_energy)) then
        recovered=tin(face)
        call gas_internal_energy(inflow_state(:,face),cell_held,accepted)
        if(.not.accepted) return
      endif
      call reactive_conserved_to_primitive(species,inflow_state(:,face),tin(face), &
        primitive,recovered,speed,accepted,internal_energy_density=cell_held)
      if(accepted) then
        value=inflow_state(:,face)
        call boundary_energy(.true.,accepted)
      end if
    end subroutine
    logical function matching_snapshot(other) result(matches)
      type(eb_reactive_tree_3d), intent(in) :: other
      integer(int64), allocatable :: other_origin(:,:),other_domain(:,:)
      integer, allocatable :: other_ratios(:)
      integer :: p
      logical :: good
      matches=.false.
      if (.not.valid_tree(other,nv)) return
      if (size(other%nodes)/=size(tree%nodes)) return
      call amr_tree_integer_layout_3d(other%topology,other_origin,other_domain,other_ratios,good)
      if (.not.good) return
      if (size(other_ratios)/=size(ratios)) return
      if (any(other_ratios/=ratios).or.any(other_domain/=domain).or.any(other_origin/=origin)) return
      do p=1,size(tree%nodes)
        if (other%topology%nodes(p)%parent/=tree%topology%nodes(p)%parent.or. &
          other%topology%nodes(p)%level/=tree%topology%nodes(p)%level) return
        if (.not.matching_eb_metrics_3d(other%nodes(p)%geometry,tree%nodes(p)%geometry)) return
      end do
      matches=.true.
    end function matching_snapshot
  end subroutine fill_reactive_eb_tree_halo_3d

  ! Provisional SDC chemistry evolves all fluid cells at each selected level,
  ! including refined parents. Reflux and restriction belong to the time driver.
  ! Estimate the first MOL predictor history without advancing public gas.
  ! Match pinned post_init ordering: root reaction, root restriction, then fine
  ! reactions. Covered solid cells never enter the reactor. Publish atomically.
  subroutine initialize_eb_tree_reaction_source_3d(species,reactions,tree,dt,subcycling,rtol,atol, &
      reaction_source,ok,evaluate_reaction,internal_energy,frame)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(eb_reactive_tree_3d), intent(in) :: tree
    real(dp), intent(in) :: dt,rtol,atol
    logical, intent(in) :: subcycling
    type(amr_tree_state_3d), intent(inout) :: reaction_source
    logical, intent(out) :: ok
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy
    type(rotating_frame_3d), optional, intent(in) :: frame
    type(amr_tree_state_3d) :: candidate,temps,rates,sources
    type(amr_tree_state_3d), allocatable :: held
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer :: node,depth,nv,n(3)
    real(dp) :: level_dt
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if(.not.valid_tree(tree,nv)) return
    if(.not.reaction_source%is_valid(tree%topology)) return
    do node=1,size(tree%nodes)
      if(size(reaction_source%fields(node)%values,1)/=nv-4) return
    enddo
    call amr_tree_integer_layout_3d(tree%topology,origin,domain,ratios,valid)
    if(.not.valid) return
    allocate(candidate%fields(size(tree%nodes)),temps%fields(size(tree%nodes)))
    do node=1,size(tree%nodes)
      n=tree%topology%nodes(node)%shape
      candidate%fields(node)%values=tree%nodes(node)%state
      allocate(temps%fields(node)%values(1,n(1),n(2),n(3)))
      temps%fields(node)%values(1,:,:,:)=tree%nodes(node)%temperature
    enddo
    rates=candidate;sources=reaction_source
    if(present(internal_energy)) then
      allocate(held);held=internal_energy
    endif
    do node=1,size(tree%nodes)
      rates%fields(node)%values=0;sources%fields(node)%values=0
    enddo
    level_dt=dt
    do depth=0,tree%topology%level_count()-1
      call react_sdc_eb_tree_3d(species,reactions,tree,rates,level_dt,rtol,atol,candidate,temps,sources,valid, &
        target_depth=depth,evaluate_reaction=evaluate_reaction,internal_energy=held, &
        old_internal_energy=internal_energy,frame=frame)
      if(.not.valid) return
      if(depth==0) then
        call average_down_reactive_eb_tree_field_3d(tree,sources,valid)
        if(.not.valid) return
      endif
      if(subcycling.and.depth<tree%topology%level_count()-1) level_dt=level_dt/real(ratios(depth+2),dp)
    enddo
    call move_alloc(sources%fields,reaction_source%fields);ok=.true.
  end subroutine initialize_eb_tree_reaction_source_3d

  subroutine react_sdc_eb_tree_3d(species,reactions,old_tree,nonreacting_rates,dt,rtol,atol, &
      state,temperature,reaction_source,ok,target_depth, &
      evaluate_reaction,frame,heat_release,internal_energy,energy_policy,energy_correction,old_internal_energy)
    type(amr_tree_state_3d), optional, intent(inout) :: heat_release,internal_energy
    type(amr_tree_state_3d), optional, intent(in) :: old_internal_energy
    type(dual_energy_policy), optional, intent(in) :: energy_policy
    real(dp), optional, intent(inout) :: energy_correction
    type(amr_tree_state_3d), allocatable :: heat,held
    type(dual_energy_policy) :: policy
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), allocatable :: potential(:,:,:)
    real(dp) :: position(3),correction,recovered,speed
    real(dp), allocatable :: node_potential(:,:,:)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(nasa7_species), intent(in) :: species(:)
    real(dp) :: primitive(reactive_nprim(size(species)))
    type(elementary_reaction), intent(in) :: reactions(:)
    type(eb_reactive_tree_3d), intent(in) :: old_tree
    type(amr_tree_state_3d), intent(in) :: nonreacting_rates
    real(dp), intent(in) :: dt,rtol,atol
    type(amr_tree_state_3d), intent(inout) :: state,temperature,reaction_source
    logical, intent(out) :: ok
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d) :: candidate,temps,sources
    integer :: nv,node,n(3),depth,i,j,k
    logical :: valid
    real(dp), allocatable :: node_heat(:,:,:),cell_held(:,:,:)
    real(dp), allocatable :: cell_old_held(:,:,:)
    ok=.false.;nv=reactive_nvar(size(species));depth=-1
    if(present(frame)) then
      if(.not.frame%is_valid()) return
      allocate(potential(1,1,1))
    end if
    if(.not.valid_tree(old_tree,nv).or.size(reactions)<1) return
    if(.not.field_shape(nonreacting_rates,nv).or..not.field_shape(reaction_source,nv-4)) return
    if(.not.all(ieee_is_finite([dt,rtol,atol])).or.min(dt,rtol,atol)<=0) return
    if(present(target_depth)) then
      depth=target_depth
      if(depth<0.or.depth>=old_tree%topology%level_count()) return
      if(.not.field_shape(state,nv).or..not.field_shape(temperature,1)) return
      candidate=state;temps=temperature
    else
      allocate(candidate%fields(size(old_tree%nodes)),temps%fields(size(old_tree%nodes)))
      do node=1,size(old_tree%nodes)
        n=old_tree%topology%nodes(node)%shape
        candidate%fields(node)%values=old_tree%nodes(node)%state
        allocate(temps%fields(node)%values(1,n(1),n(2),n(3)))
        temps%fields(node)%values(1,:,:,:)=old_tree%nodes(node)%temperature
      end do
    end if
    if(present(internal_energy)) then
      if(.not.field_shape(internal_energy,1)) return
      allocate(held);held=internal_energy
      allocate(cell_held(1,1,1))
    endif
    if(present(old_internal_energy)) then
      if(.not.present(internal_energy)) return
      if(.not.field_shape(old_internal_energy,1)) return
      allocate(cell_old_held(1,1,1))
    endif
    if(present(energy_policy).and..not.present(energy_correction)) return
    if(present(energy_correction).and..not.present(internal_energy)) return
    if(present(energy_correction)) then
      if(.not.ieee_is_finite(energy_correction)) return
      correction=energy_correction
      if(present(energy_policy)) policy=energy_policy
    endif
    if(present(heat_release)) then
      if(.not.heat_release%is_valid(old_tree%topology)) return
      if(size(heat_release%fields(1)%values,1)/=1) return
      allocate(heat);heat=heat_release
    endif
    if(present(heat_release)) allocate(node_heat(1,1,1))
    sources=reaction_source
    do node=1,size(old_tree%nodes)
      if(depth>=0.and.old_tree%topology%nodes(node)%level/=depth) cycle
      n=old_tree%topology%nodes(node)%shape
      candidate%fields(node)%values=old_tree%nodes(node)%state
      temps%fields(node)%values(1,:,:,:)=old_tree%nodes(node)%temperature
      if(present(heat_release)) heat%fields(node)%values=0
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        if(old_tree%nodes(node)%geometry%volume_fraction(i,j,k)<=0) cycle
        if(present(old_internal_energy)) cell_old_held=old_internal_energy%fields(node)%values(1,i:i,j:j,k:k)
        if(present(frame)) then
          associate(g=>old_tree%nodes(node)%geometry)
          position=[g%x_lower,g%y_lower,g%z_lower]+(real([i,j,k],dp)-.5_dp)*[g%dx,g%dy,g%dz]
          end associate
          call rotating_potential_3d(frame,position,potential(1,1,1),valid)
          if(.not.valid) return
        end if
        call react_sdc_state_3d(species,reactions,old_tree%nodes(node)%state(:,i:i,j:j,k:k), &
          old_tree%nodes(node)%temperature(i:i,j:j,k:k),nonreacting_rates%fields(node)%values(:,i:i,j:j,k:k), &
          dt,rtol,atol,candidate%fields(node)%values(:,i:i,j:j,k:k),temps%fields(node)%values(1,i:i,j:j,k:k), &
          sources%fields(node)%values(:,i:i,j:j,k:k),valid, &
          evaluate_reaction=evaluate_reaction,rotating_potential=potential,heat_release=node_heat, &
          internal_energy_density=cell_held,old_internal_energy_density=cell_old_held)
        if(.not.valid) return
        if(present(heat_release)) heat%fields(node)%values(1,i,j,k)=node_heat(1,1,1)
        if(present(internal_energy)) held%fields(node)%values(1,i,j,k)=cell_held(1,1,1)
      end do
      end do
      end do
    end do
    if(present(energy_correction)) then
      do node=1,size(old_tree%nodes)
        if(depth>=0.and.old_tree%topology%nodes(node)%level/=depth) cycle
        n=old_tree%topology%nodes(node)%shape
        if(present(frame)) then
          allocate(node_potential(n(1),n(2),n(3)),source=0._dp)
          do k=1,n(3)
          do j=1,n(2)
          do i=1,n(1)
            if(old_tree%nodes(node)%geometry%volume_fraction(i,j,k)<=0) cycle
            associate(g=>old_tree%nodes(node)%geometry)
            position=[g%x_lower,g%y_lower,g%z_lower]+(real([i,j,k],dp)-.5_dp)*[g%dx,g%dy,g%dz]
            end associate
            call rotating_potential_3d(frame,position,node_potential(i,j,k),valid)
            if(.not.valid) return
          enddo
          enddo
          enddo
        endif
        call reset_nasa7_dual_energy_patch(species,candidate%fields(node)%values, &
          held%fields(node)%values(1,:,:,:),temps%fields(node)%values(1,:,:,:), &
          old_tree%nodes(node)%geometry%volume_fraction* &
          product([old_tree%nodes(node)%geometry%dx,old_tree%nodes(node)%geometry%dy, &
          old_tree%nodes(node)%geometry%dz]),policy,correction,valid,node_potential)
        if(.not.valid) return
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(old_tree%nodes(node)%geometry%volume_fraction(i,j,k)<=0) cycle
          call recover_reactive_eb_cell_3d(species,old_tree%nodes(node)%geometry,[i,j,k], &
            candidate%fields(node)%values(:,i,j,k),temps%fields(node)%values(1,i,j,k), &
            primitive,recovered,speed,valid,frame,held%fields(node)%values(1,i,j,k))
          if(.not.valid) return
          temps%fields(node)%values(1,i,j,k)=recovered
        enddo
        enddo
        enddo
        if(allocated(node_potential)) deallocate(node_potential)
      enddo
    endif
    if(present(heat_release)) heat_release=heat
    if(present(internal_energy)) internal_energy=held
    if(present(energy_correction)) energy_correction=correction
    state=candidate;temperature=temps;reaction_source=sources;ok=.true.
  contains
    logical function field_shape(field,components) result(good)
      type(amr_tree_state_3d), intent(in) :: field
      integer, intent(in) :: components
      integer :: p
      good=.false.
      if(.not.allocated(field%fields)) return
      if(size(field%fields)/=size(old_tree%nodes)) return
      do p=1,size(field%fields)
        if(.not.allocated(field%fields(p)%values)) return
        if(any(shape(field%fields(p)%values)/=[components,old_tree%topology%nodes(p)%shape])) return
      end do
      good=.true.
    end function field_shape
  end subroutine react_sdc_eb_tree_3d

  subroutine react_reactive_eb_tree_3d(species,reactions,dt,rtol,atol,tree,chemical_change,ok,integrator)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: dt,rtol,atol
    type(eb_reactive_tree_3d), intent(inout) :: tree
    real(dp), intent(inout) :: chemical_change(:)
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: integrator
    type(eb_reactive_tree_3d) :: candidate
    logical, allocatable :: active(:,:,:)
    real(dp) :: before(size(chemical_change)),after(size(chemical_change))
    integer :: node,child,d(3),lo(3),hi(3)
    logical :: valid
    ok=.false.
    if (size(chemical_change)/=reactive_nvar(size(species))) return
    if (.not.all(ieee_is_finite([dt,rtol,atol]))) return
    if (min(dt,rtol,atol)<=0._dp.or.size(reactions)<1) return
    call integrate_reactive_eb_tree_3d(tree,before,valid)
    if (.not.valid) return
    candidate=tree
    do node=1,size(tree%nodes)
      d=tree%topology%nodes(node)%shape
      allocate(active(d(1),d(2),d(3)))
      active=tree%nodes(node)%geometry%volume_fraction>0._dp
      do child=node+1,size(tree%nodes)
        if (tree%topology%nodes(child)%parent/=node) cycle
        call bounds(tree%topology,child,lo,hi)
        active(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
      end do
      if (any(active)) then
        call advance_reactive_chemistry_3d(species,reactions,candidate%nodes(node)%state, &
          candidate%nodes(node)%temperature,d(1),d(2),d(3),dt,rtol,atol,valid,active,integrator)
        if (.not.valid) return
      end if
      deallocate(active)
    end do
    call average_down_reactive_eb_tree_3d(species,candidate,valid)
    if (.not.valid) return
    call integrate_reactive_eb_tree_3d(candidate,after,valid)
    if (.not.valid) return
    if (.not.all(ieee_is_finite(after-before))) return
    tree=candidate;chemical_change=after-before;ok=.true.
  end subroutine react_reactive_eb_tree_3d

  subroutine integrate_reactive_eb_tree_3d(tree,values,ok)
    type(eb_reactive_tree_3d), intent(in) :: tree
    real(dp), intent(inout) :: values(:)
    logical, intent(out) :: ok
    real(dp) :: candidate(size(values)),volume
    logical, allocatable :: exposed(:,:,:)
    integer :: node,child,d(3),lo(3),hi(3),i,j,k
    ok=.false.
    if (.not.valid_tree(tree,size(values))) return
    candidate=0._dp
    do node=1,size(tree%nodes)
      d=tree%topology%nodes(node)%shape
      allocate(exposed(d(1),d(2),d(3)))
      exposed=tree%nodes(node)%geometry%volume_fraction>0._dp
      do child=node+1,size(tree%nodes)
        if (tree%topology%nodes(child)%parent/=node) cycle
        call bounds(tree%topology,child,lo,hi)
        exposed(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
      end do
      volume=tree%nodes(node)%geometry%dx*tree%nodes(node)%geometry%dy*tree%nodes(node)%geometry%dz
      do k=1,d(3)
      do j=1,d(2)
      do i=1,d(1)
        if (.not.exposed(i,j,k)) cycle
        if (.not.all(ieee_is_finite(tree%nodes(node)%state(:,i,j,k)))) return
        candidate=candidate+volume*tree%nodes(node)%geometry%volume_fraction(i,j,k)*tree%nodes(node)%state(:,i,j,k)
      end do
      end do
      end do
      deallocate(exposed)
    end do
    if (.not.all(ieee_is_finite(candidate))) return
    values=candidate;ok=.true.
  end subroutine integrate_reactive_eb_tree_3d
end module eb_reactive_patch_tree_3d_mod
