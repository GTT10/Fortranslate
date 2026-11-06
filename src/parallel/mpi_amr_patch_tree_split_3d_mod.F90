module mpi_amr_patch_tree_split_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d,rotating_momentum_source_3d
  use dual_energy_mod, only: dual_energy_policy,reset_nasa7_dual_energy_cell
  use constant_volume_reactor_mod, only: forced_reactor_integrator
  use amr_patch_tree_split_3d_mod, only: amr_sdc_diffusion_evaluator
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use mpi_amr_patch_tree_chemistry_3d_mod, only: advance_mpi_amr_tree_chemistry_3d,clone_mpi_amr_tree_inputs_3d
  use mpi_amr_patch_tree_hydro_3d_mod, only: advance_mpi_amr_tree_hydro_3d
  use mpi_amr_patch_tree_transport_3d_mod, only: advance_mpi_amr_tree_transport_3d
  use mpi_amr_patch_tree_subcycle_3d_mod, only: advance_mpi_amr_tree_subcycle_hydro_3d
  use mpi_amr_patch_tree_transport_subcycle_3d_mod, only: advance_mpi_amr_tree_subcycle_transport_3d
  use reactive_3d_mod, only: compute_reactive_external_source_3d,recover_reactive_temperatures_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: recover_mpi_amr_tree_temperatures_3d
  use mpi_amr_patch_tree_reflux_3d_mod, only: reflux_mpi_amr_tree_3d
  use mpi_amr_patch_tree_hydro_3d_mod, only: compute_sdc_mpi_amr_tree_hydro_3d, compute_mol_mpi_amr_tree_hydro_3d
  use mpi_amr_patch_tree_transport_3d_mod, only: compute_sdc_mpi_amr_tree_transport_3d,agree_mpi_thermal_slip_faces
  use mpi_amr_patch_tree_chemistry_3d_mod, only: react_sdc_mpi_amr_tree_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  use mpi_amr_patch_tree_transport_3d_mod, only: agree_mpi_constant_transport_coefficients
  implicit none
  private
  public :: advance_mpi_amr_tree_split_3d,advance_mpi_amr_tree_sdc_3d
  public :: advance_mpi_amr_tree_sdc_subcycle_3d
  public :: agree_mpi_amr_tree_energy_policy_3d,close_mpi_amr_tree_held_energy_3d
  interface agree_mpi_amr_tree_energy_policy_3d
    module procedure agree_tree_energy_policy
  end interface
  abstract interface
    subroutine mpi_amr_sdc_auxiliary_source(state,temperature,interval,stage,depth,rate,comm,ok)
      import :: amr_tree_state_3d,dp,MPI_Comm
      type(amr_tree_state_3d), intent(in) :: state,temperature
      real(dp), intent(in) :: interval
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
    end subroutine
    ! Time denotes the level interval start relative to the outer step.
    subroutine mpi_amr_sdc_recursive_source(state,temperature,interval,stage,depth,time,rate,comm,ok)
      import :: amr_tree_state_3d,dp,MPI_Comm
      type(amr_tree_state_3d), intent(in) :: state,temperature
      real(dp), intent(in) :: interval,time
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
    end subroutine
    subroutine mpi_amr_sdc_recursive_halo_source(state,temperature,interval,stage,depth,time,rate,comm,ok, &
        ancestor_start,ancestor_end,ancestor_fraction)
      import :: amr_tree_state_3d,dp,MPI_Comm
      type(amr_tree_state_3d), intent(in) :: state,temperature,ancestor_start,ancestor_end
      real(dp), intent(in) :: interval,time,ancestor_fraction(:)
      integer, intent(in) :: stage,depth
      type(amr_tree_state_3d), intent(inout) :: rate
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
    end subroutine
    subroutine mpi_amr_sdc_post_level(depth,time,interval,iteration,ncycle,comm,ok)
      import :: dp,MPI_Comm
      integer, intent(in) :: depth,iteration,ncycle
      real(dp), intent(in) :: time,interval
      type(MPI_Comm), intent(in) :: comm
      logical, intent(out) :: ok
    end subroutine
  end interface
  public :: mpi_amr_sdc_auxiliary_source,mpi_amr_sdc_recursive_source
  public :: mpi_amr_sdc_recursive_halo_source,mpi_amr_sdc_post_level
contains
  subroutine advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,state,temperature,spacing,dt,faces, &
      iterations,chemistry,transport_flags,force,rtol,atol,boundary_change,source_change,maximum_correction,minimum_theta, &
      comm,ok,reaction_source,inflow_state,inflow_temperature,wall_temperature,wall_velocity,outlet_pressure, &
      inflow_mach,inflow_pressure,difmag,level_steps,evaluate_diffusion,evaluate_reaction, &
      constant_coefficients,auxiliary_source,auxiliary_change,halo_source,post_level,heat_release, &
      mol_mode,mol_order,mol_lax_friedrichs, &
      subcycling,integrated_flux,old_predictor,new_predictor,internal_energy,energy_correction,energy_policy, &
      internal_energy_integrated_flux,physical_integrated_flux,physical_internal_energy_integrated_flux,frame,domain_lower, &
      auxiliary_samples_internal_energy,internal_energy_pressure_work_change,isothermal_slip_faces)
    logical, optional, intent(in) :: isothermal_slip_faces(6)
    logical :: thermal_slip_check(6)
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3)
    logical, optional, intent(in) :: subcycling
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: integrated_flux(:)
    type(amr_tree_state_3d), optional, intent(inout) :: old_predictor,new_predictor
    logical :: cycle_levels,mol
    ! Opt-in callback input: physical gas followed by actual rho e. Rates
    ! retain the physical gas width; particles/films have zero direct UEINT source.
    logical, optional, intent(in) :: auxiliary_samples_internal_energy
    logical :: sample_held
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy
    ! Signed pressure-work increment; commit only with the outer gas/H transaction.
    real(dp), optional, intent(inout) :: internal_energy_pressure_work_change
    real(dp), optional, intent(inout) :: energy_correction
    type(dual_energy_policy), optional, intent(in) :: energy_policy
    type(dual_energy_policy) :: policy
    logical :: stored
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: internal_energy_integrated_flux(:)
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: physical_integrated_flux(:)
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: physical_internal_energy_integrated_flux(:)
    logical, optional, intent(in) :: mol_mode,mol_lax_friedrichs
    integer, optional, intent(in) :: mol_order
    type(amr_tree_state_3d), optional, intent(inout) :: heat_release
    type(amr_tree_state_3d), allocatable :: heat
    procedure(mpi_amr_sdc_recursive_source), optional :: auxiliary_source
    procedure(mpi_amr_sdc_recursive_halo_source), optional :: halo_source
    procedure(mpi_amr_sdc_post_level), optional :: post_level
    real(dp), optional, intent(inout) :: auxiliary_change(:)
    real(dp), optional, intent(in) :: constant_coefficients(:)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    procedure(amr_sdc_diffusion_evaluator), optional :: evaluate_diffusion
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: spacing(3),dt,force(3),rtol,atol
    character(len=*), intent(in) :: faces(6)
    integer, intent(in) :: iterations
    logical, intent(in) :: chemistry,transport_flags(4)
    real(dp), intent(inout) :: boundary_change(:),source_change(:),maximum_correction,minimum_theta
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6),difmag
    integer(int64), optional, intent(inout) :: level_steps(:)
    type(amr_tree_state_3d) :: candidate,temps,old_endpoint,new_endpoint,source0,source1
    type(amr_tree_state_3d), allocatable :: ir,held
    type(amr_tree_flux_3d), allocatable :: integrated(:),integrated_aux(:),integrated_physical(:),integrated_aux_physical(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:),counts(:),executed(:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: starts(:),finishes(:),frame_h(:,:)
    real(dp) :: budget(size(boundary_change)),volume(size(source_change)),largest,theta
    real(dp) :: auxiliary(size(source_change)),active_time,energy_delta,reset_change
    real(dp) :: pressure_delta
    real(dp), allocatable :: step_auxiliary(:)
    integer :: nv,node,n(3),level,rank,status,flags(19),low(19),high(19)
    logical :: valid,has_auxiliary,preserve_mol_old
    integer :: heat_supplied,heat_min,heat_max,heat_status
    has_auxiliary=present(auxiliary_source).or.present(halo_source)
    ok=.false.;nv=5+size(species)
    call agree_mpi_thermal_slip_faces(isothermal_slip_faces,comm,thermal_slip_check,ok, &
      faces,transport_flags(2),wall_temperature)
    if(.not.ok) return
    ok=.false.
    mol=.false.
    if(present(mol_mode)) mol=mol_mode
    cycle_levels=.true.
    if(present(subcycling)) cycle_levels=subcycling
    sample_held=.false.
    if(present(auxiliary_samples_internal_energy)) sample_held=auxiliary_samples_internal_energy
    call agree_tree_energy_policy(mol,internal_energy,energy_correction,energy_policy,comm,policy,stored,valid)
    if(.not.valid) return
    heat_supplied=merge(1,0,present(heat_release))
    call MPI_Allreduce(heat_supplied,heat_min,1,MPI_INTEGER,MPI_MIN,comm,heat_status)
    if(heat_status/=MPI_SUCCESS) return
    call MPI_Allreduce(heat_supplied,heat_max,1,MPI_INTEGER,MPI_MAX,comm,heat_status)
    if(heat_status/=MPI_SUCCESS.or.heat_min/=heat_max) return
    call agree_mpi_constant_transport_coefficients(constant_coefficients,comm,ok,species=species)
    if(.not.ok) return
    ok=.false.
    flags=[merge(1,0,chemistry),merge(1,0,present(reaction_source)),merge(1,0,present(level_steps)), &
      merge(1,0,present(evaluate_reaction)),merge(1,0,present(auxiliary_source)), &
      merge(1,0,present(halo_source)),merge(1,0,present(post_level)),merge(1,0,cycle_levels), &
      merge(1,0,present(integrated_flux)),merge(1,0,present(old_predictor)),merge(1,0,present(new_predictor)), &
      merge(1,0,present(internal_energy_integrated_flux)),merge(1,0,present(physical_integrated_flux)), &
      merge(1,0,present(physical_internal_energy_integrated_flux)), &
      merge(1,0,present(frame)),merge(1,0,present(domain_lower)), &
      merge(1,0,present(auxiliary_samples_internal_energy)),merge(1,0,sample_held), &
      merge(1,0,present(internal_energy_pressure_work_change))]
    call MPI_Allreduce(flags,low,size(flags),MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,size(flags),MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(low/=high)) return
    if(present(internal_energy_pressure_work_change).and..not.stored) return
    valid=present(frame).eqv.present(domain_lower)
    call agree(valid)
    if(.not.valid) return
    if(present(frame)) then
      valid=mol.and.present(domain_lower).and.frame%is_valid()
      if(present(domain_lower)) valid=valid.and.all(ieee_is_finite(domain_lower))
      call agree(valid)
      if(.not.valid) return
    endif
    valid=size(boundary_change)==nv.and.size(source_change)==nv
    valid=valid.and.all(ieee_is_finite([spacing,dt])).and.minval(spacing)>0.and.dt>0
    valid=valid.and.(chemistry.or..not.present(reaction_source))
    valid=valid.and.(chemistry.or..not.present(evaluate_reaction))
    valid=valid.and.(.not.sample_held.or.(stored.and.has_auxiliary))
    valid=valid.and.(has_auxiliary.eqv.present(auxiliary_change))
    valid=valid.and..not.(present(auxiliary_source).and.present(halo_source))
    if(present(auxiliary_change)) valid=valid.and.size(auxiliary_change)==nv
    call agree(valid)
    if(.not.valid) return
    call clone_mpi_amr_tree_inputs_3d(tree,owners,size(species),state,temperature,candidate,temps,comm,rank,valid)
    if(.not.valid) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid)
    if(.not.valid) return
    if(present(frame)) then
      allocate(frame_h(3,size(tree%nodes)))
      frame_h(:,1)=spacing
      do node=2,size(tree%nodes)
        frame_h(:,node)=frame_h(:,tree%nodes(node)%parent)/real(tree%nodes(node)%patch%refinement_ratio,dp)
      enddo
    endif
    allocate(counts(size(ratios)),executed(size(ratios)),starts(size(ratios)),finishes(size(ratios)))
    counts(1)=1
    do level=2,size(ratios)
      if(counts(level-1)>huge(1_int64)/int(ratios(level),int64)) return
      counts(level)=counts(level-1)
      if(cycle_levels) counts(level)=counts(level)*int(ratios(level),int64)
    end do
    if(present(level_steps)) then
      valid=size(level_steps)==size(counts)
      call agree(valid)
      if(.not.valid) return
    end if
    if(present(heat_release)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,heat_release,1,comm,valid)
      if(.not.valid) return
      allocate(heat);heat=heat_release
    endif
    if(present(reaction_source)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,reaction_source,nv-4,comm,valid)
      if(.not.valid) return
    end if
    preserve_mol_old=.false.
    if(present(mol_mode)) preserve_mol_old=mol_mode.and.stored.and..not.present(frame)
    preserve_mol_old=preserve_mol_old.and.policy%update_total_from_internal
    energy_delta=0;pressure_delta=0
    if(stored) then
      call validate_mpi_amr_tree_local_3d(tree,owners,internal_energy,1,comm,valid)
      if(.not.valid) return
      allocate(held);held=internal_energy
      if(.not.preserve_mol_old) then
        call average_down_mpi_amr_tree_3d(tree,owners,held,1,comm,valid)
        if(.not.valid) return
      endif
    endif
    if(preserve_mol_old) then
      call validate_mpi_mol_old_state(species,tree,owners,candidate,temps,held,comm,valid)
    else
      call average_down_mpi_amr_tree_3d(tree,owners,candidate,nv,comm,valid)
      if(.not.valid) return
      call recover_from(0,valid)
    endif
    if(.not.valid) return
    allocate(old_endpoint%fields(size(tree%nodes)),new_endpoint%fields(size(tree%nodes)))
    source0=state;source1=state
    if(chemistry) allocate(ir)
    if(chemistry) allocate(ir%fields(size(tree%nodes)))
    do node=1,size(tree%nodes)
      if(owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      allocate(old_endpoint%fields(node)%values(nv+1+merge(1,0,stored),n(1),n(2),n(3)))
      old_endpoint%fields(node)%values=0
      new_endpoint%fields(node)%values=old_endpoint%fields(node)%values
      source0%fields(node)%values=0;source1%fields(node)%values=0
      if(chemistry) then
        allocate(ir%fields(node)%values(nv-4,n(1),n(2),n(3)));ir%fields(node)%values=0
      end if
    end do
    if(present(reaction_source)) ir=reaction_source
    if(has_auxiliary) then
      allocate(step_auxiliary(nv))
    endif
    auxiliary=0;active_time=0
    starts=0;finishes=1;executed=0;budget=0;volume=0;largest=0;theta=1
    call interval(0,0._dp,1._dp,1,1,integrated,valid,integrated_aux,integrated_physical,integrated_aux_physical)
    if(.not.valid.or.any(executed/=counts)) return
    if(.not.all(ieee_is_finite(budget)).or..not.all(ieee_is_finite(volume))) return
    if(.not.all(ieee_is_finite(auxiliary))) return
    if(present(mol_mode)) then
      if(mol_mode) then
        valid=.true.
        do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
          if(any(candidate%fields(node)%values(6:,:,:,:)<0._dp)) valid=.false.
        enddo
        call agree(valid)
        if(.not.valid) return
      endif
    endif
    if(stored) then
      valid=ieee_is_finite(energy_correction+energy_delta);call agree(valid)
      if(.not.valid) return
      internal_energy=held;energy_correction=energy_correction+energy_delta
      if(present(internal_energy_integrated_flux)) internal_energy_integrated_flux=integrated_aux
      if(present(physical_internal_energy_integrated_flux)) &
        physical_internal_energy_integrated_flux=integrated_aux_physical
    endif
    if(present(internal_energy_pressure_work_change)) internal_energy_pressure_work_change=pressure_delta
    state=candidate;temperature=temps;boundary_change=budget;source_change=volume
    if(present(auxiliary_change)) auxiliary_change=auxiliary
    maximum_correction=largest;minimum_theta=theta
    if(present(heat_release)) heat_release=heat
    if(present(reaction_source)) reaction_source=ir
    if(present(level_steps)) level_steps=executed
    if(present(integrated_flux)) integrated_flux=integrated
    if(present(physical_integrated_flux)) physical_integrated_flux=integrated_physical
    if(present(old_predictor)) old_predictor=source0
    if(present(new_predictor)) new_predictor=source1
    ok=.true.
  contains
    subroutine source_adapter(stage_state,stage_temperature,interval_dt,stage,depth,rate,source_comm,success)
      type(amr_tree_state_3d), intent(in) :: stage_state,stage_temperature
      type(amr_tree_state_3d), intent(inout) :: rate
      real(dp), intent(in) :: interval_dt
      integer, intent(in) :: stage,depth
      type(MPI_Comm), intent(in) :: source_comm
      logical, intent(out) :: success
      real(dp) :: fractions(depth),source_time
      integer :: ancestor
      source_time=active_time
      if(stage>0) source_time=source_time+interval_dt
      if(present(halo_source)) then
        do ancestor=1,depth
          fractions(ancestor)=max(0._dp,min(1._dp, &
            (source_time/dt-starts(ancestor))/(finishes(ancestor)-starts(ancestor))))
        enddo
        call halo_source(stage_state,stage_temperature,interval_dt,stage,depth,active_time,rate,source_comm,success, &
          old_endpoint,new_endpoint,fractions)
      else
        call auxiliary_source(stage_state,stage_temperature,interval_dt,stage,depth,active_time,rate,source_comm,success)
      endif
    end subroutine source_adapter
    recursive subroutine interval(depth,begin_time,end_time,iteration,ncycle,effective,success, &
        auxiliary_effective,physical_effective,auxiliary_physical_effective)
      integer, intent(in) :: depth,iteration,ncycle
      real(dp), intent(in) :: begin_time,end_time
      type(amr_tree_flux_3d), allocatable, intent(out) :: effective(:)
      type(amr_tree_flux_3d), allocatable, optional, intent(out) :: auxiliary_effective(:)
      type(amr_tree_flux_3d), allocatable, intent(out) :: physical_effective(:)
      type(amr_tree_flux_3d), allocatable, intent(out) :: auxiliary_physical_effective(:)
      logical, intent(out) :: success
      type(amr_tree_flux_3d), allocatable :: total(:),child_flux(:),total_aux(:),child_aux(:),effective_aux(:), &
        total_physical(:),child_physical(:),total_aux_physical(:),child_aux_physical(:)
      type(amr_tree_state_3d) :: op,np,registration_state,registered_held
      real(dp), allocatable :: local_energy_delta
      type(dual_energy_policy), allocatable :: local_policy
      real(dp), allocatable :: local_pressure_delta
      real(dp) :: fraction0(depth),fraction1(depth),change(nv),external(nv),correction(nv),registered_correction(nv)
      real(dp) :: max_change,local_theta,sub_begin,sub_end,local_dt,scalar_correction(1),scalar_max,registered_max
      integer :: peer,ancestor,substep,ratio
      logical :: good
      success=.false.;local_dt=dt*(end_time-begin_time)
      if(present(internal_energy_pressure_work_change)) allocate(local_pressure_delta,source=0._dp)
      if(stored) then
        allocate(local_energy_delta,source=0._dp);allocate(local_policy,source=policy)
      endif
      if(end_time<=begin_time.or.local_dt<=0) return
      starts(depth+1)=begin_time;finishes(depth+1)=end_time
      do ancestor=1,depth
        fraction0(ancestor)=max(0._dp,min(1._dp,(begin_time-starts(ancestor))/(finishes(ancestor)-starts(ancestor))))
        fraction1(ancestor)=max(0._dp,min(1._dp,(end_time-starts(ancestor))/(finishes(ancestor)-starts(ancestor))))
      end do
      active_time=dt*begin_time
      call pack_depth(depth,old_endpoint)
      if(depth==0) then
        if(has_auxiliary) then
        call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,candidate,temps,spacing,local_dt,faces, &
          iterations,chemistry,transport_flags,force,rtol,atol,change,external,max_change,local_theta,comm,good, &
          ir,inflow_state,inflow_temperature,wall_temperature,wall_velocity,outlet_pressure,inflow_mach,inflow_pressure, &
          difmag,target_depth=depth,integrated_flux=total,old_predictor=op,new_predictor=np, &
          evaluate_diffusion=evaluate_diffusion,evaluate_reaction=evaluate_reaction, &
          constant_coefficients=constant_coefficients,auxiliary_source=source_adapter,auxiliary_change=step_auxiliary, &
          heat_release=heat, &
          mol_mode=mol_mode,mol_order=mol_order,mol_lax_friedrichs=mol_lax_friedrichs, &
          internal_energy=held,energy_correction=local_energy_delta,energy_policy=local_policy, &
          internal_energy_integrated_flux=total_aux,physical_integrated_flux=total_physical, &
          physical_internal_energy_integrated_flux=total_aux_physical,frame=frame,domain_lower=domain_lower, &
          internal_energy_pressure_work_change=local_pressure_delta, &
          auxiliary_samples_internal_energy=sample_held, &
          isothermal_slip_faces=isothermal_slip_faces)
        else
        call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,candidate,temps,spacing,local_dt,faces, &
          iterations,chemistry,transport_flags,force,rtol,atol,change,external,max_change,local_theta,comm,good, &
          ir,inflow_state,inflow_temperature,wall_temperature,wall_velocity,outlet_pressure,inflow_mach,inflow_pressure, &
          difmag,target_depth=depth,integrated_flux=total,old_predictor=op,new_predictor=np, &
          evaluate_diffusion=evaluate_diffusion,evaluate_reaction=evaluate_reaction, &
          constant_coefficients=constant_coefficients,heat_release=heat, &
          mol_mode=mol_mode,mol_order=mol_order,mol_lax_friedrichs=mol_lax_friedrichs, &
          internal_energy=held,energy_correction=local_energy_delta,energy_policy=local_policy, &
          internal_energy_integrated_flux=total_aux,physical_integrated_flux=total_physical, &
          physical_internal_energy_integrated_flux=total_aux_physical,frame=frame,domain_lower=domain_lower, &
          internal_energy_pressure_work_change=local_pressure_delta, &
          isothermal_slip_faces=isothermal_slip_faces)
        endif
      else
        if(has_auxiliary) then
        call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,candidate,temps,spacing,local_dt,faces, &
          iterations,chemistry,transport_flags,force,rtol,atol,change,external,max_change,local_theta,comm,good, &
          ir,inflow_state,inflow_temperature,wall_temperature,wall_velocity,outlet_pressure,inflow_mach,inflow_pressure, &
          difmag,target_depth=depth,coarse_start=old_endpoint,coarse_end=new_endpoint, &
          fractions_start=fraction0,fractions_end=fraction1,source_start=source0,source_end=source1, &
          integrated_flux=total,old_predictor=op,new_predictor=np, &
          evaluate_diffusion=evaluate_diffusion,evaluate_reaction=evaluate_reaction, &
          constant_coefficients=constant_coefficients,auxiliary_source=source_adapter,auxiliary_change=step_auxiliary, &
          heat_release=heat, &
          mol_mode=mol_mode,mol_order=mol_order,mol_lax_friedrichs=mol_lax_friedrichs, &
          internal_energy=held,energy_correction=local_energy_delta,energy_policy=local_policy, &
          internal_energy_integrated_flux=total_aux,physical_integrated_flux=total_physical, &
          physical_internal_energy_integrated_flux=total_aux_physical,frame=frame,domain_lower=domain_lower, &
          internal_energy_pressure_work_change=local_pressure_delta, &
          auxiliary_samples_internal_energy=sample_held, &
          isothermal_slip_faces=isothermal_slip_faces)
        else
        call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,candidate,temps,spacing,local_dt,faces, &
          iterations,chemistry,transport_flags,force,rtol,atol,change,external,max_change,local_theta,comm,good, &
          ir,inflow_state,inflow_temperature,wall_temperature,wall_velocity,outlet_pressure,inflow_mach,inflow_pressure, &
          difmag,target_depth=depth,coarse_start=old_endpoint,coarse_end=new_endpoint, &
          fractions_start=fraction0,fractions_end=fraction1,source_start=source0,source_end=source1, &
          integrated_flux=total,old_predictor=op,new_predictor=np, &
          evaluate_diffusion=evaluate_diffusion,evaluate_reaction=evaluate_reaction, &
          constant_coefficients=constant_coefficients,heat_release=heat, &
          mol_mode=mol_mode,mol_order=mol_order,mol_lax_friedrichs=mol_lax_friedrichs, &
          internal_energy=held,energy_correction=local_energy_delta,energy_policy=local_policy, &
          internal_energy_integrated_flux=total_aux,physical_integrated_flux=total_physical, &
          physical_internal_energy_integrated_flux=total_aux_physical,frame=frame,domain_lower=domain_lower, &
          internal_energy_pressure_work_change=local_pressure_delta, &
          isothermal_slip_faces=isothermal_slip_faces)
        endif
      end if
      if(.not.good) return
      if(stored) energy_delta=energy_delta+local_energy_delta
      if(present(internal_energy_pressure_work_change)) then
        if(.not.ieee_is_finite(local_pressure_delta)) return
        if(local_pressure_delta>0) then
          if(pressure_delta>huge(1._dp)-local_pressure_delta) return
        else if(local_pressure_delta<0) then
          if(pressure_delta< -huge(1._dp)-local_pressure_delta) return
        endif
        pressure_delta=pressure_delta+local_pressure_delta
      endif
      if(has_auxiliary) auxiliary=auxiliary+step_auxiliary
      executed(depth+1)=executed(depth+1)+1
      theta=min(theta,local_theta);volume=volume+external
      if(depth==0) budget=budget+change
      call pack_depth(depth,new_endpoint)
      do peer=1,size(tree%nodes)
        if(owners(peer)/=rank) cycle
        if(tree%nodes(peer)%level/=depth) cycle
        source0%fields(peer)%values=op%fields(peer)%values
        source1%fields(peer)%values=np%fields(peer)%values
      end do
      if(depth+1<size(counts)) then
        ratio=1
        if(cycle_levels) ratio=ratios(depth+2)
        sub_begin=begin_time
        do substep=1,ratio
          sub_end=begin_time+(end_time-begin_time)*real(substep,dp)/real(ratio,dp)
          if(substep==ratio) sub_end=end_time
          call interval(depth+1,sub_begin,sub_end,substep,ratio,child_flux,good,child_aux,child_physical,child_aux_physical)
          if(.not.good) return
          do peer=1,size(tree%nodes)
            if(owners(peer)/=rank) cycle
            if(tree%nodes(peer)%level<=depth) cycle
            if(.not.allocated(total(peer)%x)) then
              total(peer)=child_flux(peer)
              total_physical(peer)=child_physical(peer)
            else
              total(peer)%x=total(peer)%x+child_flux(peer)%x
              total(peer)%y=total(peer)%y+child_flux(peer)%y
              total(peer)%z=total(peer)%z+child_flux(peer)%z
              total_physical(peer)%x=total_physical(peer)%x+child_physical(peer)%x
              total_physical(peer)%y=total_physical(peer)%y+child_physical(peer)%y
              total_physical(peer)%z=total_physical(peer)%z+child_physical(peer)%z
            end if
          end do
          if(stored) then
            do peer=1,size(tree%nodes)
              if(owners(peer)/=rank) cycle
              if(tree%nodes(peer)%level<=depth) cycle
              if(.not.allocated(total_aux(peer)%x)) then
                total_aux(peer)=child_aux(peer)
                total_aux_physical(peer)=child_aux_physical(peer)
              else
                total_aux(peer)%x=total_aux(peer)%x+child_aux(peer)%x
                total_aux(peer)%y=total_aux(peer)%y+child_aux(peer)%y
                total_aux(peer)%z=total_aux(peer)%z+child_aux(peer)%z
                total_aux_physical(peer)%x=total_aux_physical(peer)%x+child_aux_physical(peer)%x
                total_aux_physical(peer)%y=total_aux_physical(peer)%y+child_aux_physical(peer)%y
                total_aux_physical(peer)%z=total_aux_physical(peer)%z+child_aux_physical(peer)%z
              endif
            enddo
          endif
          sub_begin=sub_end
        end do
      end if
      if(mol.and..not.chemistry.and.iterations>1) then
        ! Reflux the accepted two-half-step state before canonicalizing the
        ! pinned old-half register on a private copy for the parent interval.
        if(stored) then
        call reflux_mpi_amr_tree_3d(tree,owners,candidate,nv,total_physical,spacing,faces,correction,max_change,comm,good, &
          minimum_depth=depth,synchronized_flux=physical_effective,first_species=6,last_species=nv, &
          auxiliary_state=held,auxiliary_flux=total_aux_physical,auxiliary_synchronized_flux=auxiliary_physical_effective)
        else
        call reflux_mpi_amr_tree_3d(tree,owners,candidate,nv,total_physical,spacing,faces,correction,max_change,comm,good, &
          minimum_depth=depth,synchronized_flux=physical_effective,first_species=merge(6,0,mol), &
          last_species=merge(nv,0,mol))
        endif
        if(.not.good) return
        registration_state=candidate
        if(stored) then
        registered_held=held
        call reflux_mpi_amr_tree_3d(tree,owners,registration_state,nv,total,spacing,faces,registered_correction, &
          registered_max,comm,good,minimum_depth=depth,synchronized_flux=effective,auxiliary_state=registered_held, &
          auxiliary_flux=total_aux,auxiliary_synchronized_flux=effective_aux)
        else
        call reflux_mpi_amr_tree_3d(tree,owners,registration_state,nv,total,spacing,faces,registered_correction, &
          registered_max,comm,good,minimum_depth=depth,synchronized_flux=effective)
        endif
      else
        if(stored) then
        call reflux_mpi_amr_tree_3d(tree,owners,candidate,nv,total,spacing,faces,correction,max_change,comm,good, &
          minimum_depth=depth,synchronized_flux=effective,first_species=6,last_species=nv, &
          auxiliary_state=held,auxiliary_flux=total_aux,auxiliary_synchronized_flux=effective_aux)
        if(good) auxiliary_physical_effective=effective_aux
        else
        call reflux_mpi_amr_tree_3d(tree,owners,candidate,nv,total,spacing,faces,correction,max_change,comm,good, &
          minimum_depth=depth,synchronized_flux=effective,first_species=merge(6,0,mol),last_species=merge(nv,0,mol))
        endif
        if(good) physical_effective=effective
      endif
      if(.not.good) return
      largest=max(largest,max_change)
      if(depth==0) budget=budget+correction
      if(stored) then
        if(present(auxiliary_effective)) auxiliary_effective=effective_aux
      endif
      call recover_from(depth,good)
      if(.not.good) return
      if(chemistry.and.depth+1<tree%level_count()) then
        call average_down_mpi_amr_tree_3d(tree,owners,ir,nv-4,comm,good,minimum_depth=depth)
        if(.not.good) return
      endif
      if(present(heat_release).and.depth+1<tree%level_count()) then
        call average_down_mpi_amr_tree_3d(tree,owners,heat,1,comm,good,minimum_depth=depth)
        if(.not.good) return
      endif
      if(present(post_level)) then
        call post_level(depth,dt*begin_time,local_dt,iteration,ncycle,comm,good)
        call agree(good)
        if(.not.good) return
      endif
      success=.true.
    end subroutine interval
    subroutine pack_depth(depth,packed)
      integer, intent(in) :: depth
      type(amr_tree_state_3d), intent(inout) :: packed
      integer :: peer
      do peer=1,size(tree%nodes)
        if(owners(peer)/=rank) cycle
        if(tree%nodes(peer)%level/=depth) cycle
        packed%fields(peer)%values(:nv,:,:,:)=candidate%fields(peer)%values
        packed%fields(peer)%values(nv+1,:,:,:)=temps%fields(peer)%values(1,:,:,:)
        if(stored) packed%fields(peer)%values(nv+2,:,:,:)=held%fields(peer)%values(1,:,:,:)
      end do
    end subroutine pack_depth
    subroutine recover_from(depth,success)
      integer, intent(in) :: depth
      logical, intent(out) :: success
      integer :: selected,peer,dims(3),i,j,k
      real(dp), allocatable :: recovered(:,:,:),potential(:,:,:)
      real(dp) :: position(3)
      if(stored) then
        call reset_mpi_tree_energy(species,tree,owners,candidate,temps,held,spacing,policy,reset_change,comm,success, &
          minimum_depth=depth,frame=frame,domain_lower=domain_lower)
        if(.not.success) return
        if(present(frame)) then
          call average_down_mpi_amr_tree_3d(tree,owners,candidate,nv,comm,success,minimum_depth=depth)
          if(.not.success) return
          call average_down_mpi_amr_tree_3d(tree,owners,held,1,comm,success,minimum_depth=depth)
          if(.not.success) return
          call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,candidate,temps,comm,success, &
            internal_energy=held)
          if(.not.success) return
        endif
        energy_delta=energy_delta+reset_change
        return
      endif
      if(present(frame)) then
        success=.true.
        do peer=1,size(tree%nodes)
          if(owners(peer)/=rank.or.tree%nodes(peer)%level<depth) cycle
          dims=tree%nodes(peer)%shape
          allocate(recovered(dims(1),dims(2),dims(3)),potential(dims(1),dims(2),dims(3)))
          do k=1,dims(3)
          do j=1,dims(2)
          do i=1,dims(1)
            position=domain_lower+(real(origin(:,peer),dp)+real([i,j,k],dp)-.5_dp)*frame_h(:,peer)
            call rotating_potential_3d(frame,position,potential(i,j,k),success)
            if(.not.success) exit
          enddo
          if(.not.success) exit
          enddo
          if(.not.success) exit
          enddo
          if(success) then
            call recover_reactive_temperatures_3d(species,candidate%fields(peer)%values,temps%fields(peer)%values(1,:,:,:), &
              dims(1),dims(2),dims(3),recovered,success,rotating_potential=potential)
          endif
          if(.not.success) exit
          temps%fields(peer)%values(1,:,:,:)=recovered
          deallocate(recovered,potential)
        enddo
        call agree(success)
        return
      endif
      do selected=depth,tree%level_count()-1
        call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,candidate,temps,comm,success,selected)
        if(.not.success) return
      end do
      success=.true.
    end subroutine recover_from
    subroutine agree(condition)
      logical, intent(inout) :: condition
      logical :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      condition=.false.
      if(ierr==MPI_SUCCESS) condition=result
    end subroutine agree
  end subroutine advance_mpi_amr_tree_sdc_subcycle_3d

  ! Whole-tree MOL uses ordered level advances with equal dt; CTU retains common stages.
  recursive subroutine advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,state,temperature,spacing,dt,faces, &
      iterations,chemistry,transport_flags,force,rtol,atol,boundary_change,source_change,maximum_correction,minimum_theta, &
      comm,ok,reaction_source,inflow_state,inflow_temperature,wall_temperature,wall_velocity,outlet_pressure, &
      inflow_mach,inflow_pressure,difmag,target_depth,coarse_start,coarse_end,fractions_start,fractions_end, &
      source_start,source_end,integrated_flux,old_predictor,new_predictor,evaluate_diffusion,evaluate_reaction, &
      constant_coefficients,auxiliary_source,auxiliary_change,heat_release,mol_mode,mol_order,mol_lax_friedrichs, &
      internal_energy,energy_correction,energy_policy,internal_energy_integrated_flux,physical_integrated_flux, &
      physical_internal_energy_integrated_flux,frame,domain_lower, &
      auxiliary_samples_internal_energy,internal_energy_pressure_work_change,isothermal_slip_faces)
    logical, optional, intent(in) :: isothermal_slip_faces(6)
    logical :: thermal_slip_check(6)
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3)
    ! Opt-in callback input: physical gas followed by actual rho e. Rates
    ! retain the physical gas width; particles/films have zero direct UEINT source.
    logical, optional, intent(in) :: auxiliary_samples_internal_energy
    logical :: sample_held
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy
    ! Signed pressure-work increment; commit only with the outer gas/H transaction.
    real(dp), optional, intent(inout) :: internal_energy_pressure_work_change
    real(dp), optional, intent(inout) :: energy_correction
    type(dual_energy_policy), optional, intent(in) :: energy_policy
    type(dual_energy_policy) :: policy
    logical :: stored
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: internal_energy_integrated_flux(:)
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: physical_integrated_flux(:)
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: physical_internal_energy_integrated_flux(:)
    logical, optional, intent(in) :: mol_mode,mol_lax_friedrichs
    integer, optional, intent(in) :: mol_order
    type(amr_tree_state_3d), optional, intent(inout) :: heat_release
    type(amr_tree_state_3d), allocatable :: heat
    real(dp), optional, intent(in) :: constant_coefficients(:)
    procedure(mpi_amr_sdc_auxiliary_source), optional :: auxiliary_source
    real(dp), optional, intent(inout) :: auxiliary_change(:)
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    procedure(amr_sdc_diffusion_evaluator), optional :: evaluate_diffusion
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: spacing(3),dt,force(3),rtol,atol
    character(len=*), intent(in) :: faces(6)
    integer, intent(in) :: iterations
    logical, intent(in) :: chemistry,transport_flags(4)
    real(dp), intent(inout) :: boundary_change(:),source_change(:),maximum_correction,minimum_theta
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6),difmag
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d), optional, intent(in) :: coarse_start,coarse_end,source_start,source_end
    real(dp), optional, intent(in) :: fractions_start(:),fractions_end(:)
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: integrated_flux(:)
    type(amr_tree_state_3d), optional, intent(inout) :: old_predictor,new_predictor
    type(amr_tree_state_3d) :: mol_old_hydro,mol_new_hydro,old_arate,new_arate,predictor_gas,predictor_held,predictor_reset
    type(amr_tree_state_3d), allocatable :: old_held,held,transport_start,transport_end
    type(amr_tree_flux_3d), allocatable :: old_aflux(:),new_aflux(:),aflux(:),physical_aflux(:)
    real(dp) :: initial_energy_delta,predictor_energy_delta,final_energy_delta,energy_delta
    real(dp) :: pressure_delta,pressure_integral(1)
    type(amr_tree_state_3d), allocatable :: pressure_work
    type(amr_tree_flux_3d), allocatable :: mol_old_flux(:),mol_new_flux(:),registration(:)
    type(amr_tree_flux_3d), allocatable :: effective_mol_flux(:),effective_mol_aux(:)
    integer(int64), allocatable :: frame_origin(:,:),frame_domain(:,:)
    integer, allocatable :: frame_ratios(:)
    real(dp), allocatable :: frame_h(:,:)
    integer :: reconstruction_order
    logical :: mol,lf
    real(dp), allocatable :: midpoint(:)
    integer :: depth,selection(14),selection_low(14),selection_high(14)
    type(amr_tree_state_3d) :: op,np,old_aux,new_aux
    type(amr_tree_state_3d) :: old,ot,candidate,ct,old_ext,new_ext,old_tr,new_tr,hydro,rates,ir
    type(amr_tree_flux_3d), allocatable :: flux(:),old_flux(:),new_flux(:)
    real(dp) :: budget(size(boundary_change)),volume(size(source_change)),correction(size(boundary_change))
    real(dp) :: first(size(source_change)),last(size(source_change)),largest,theta,stage_theta
    real(dp) :: auxiliary(size(source_change))
    integer :: nv,node,n(3),iteration,c,rank,status,ios
    integer :: switches(19),low_switches(19),high_switches(19)
    real(dp) :: controls(9),low_controls(9),high_controls(9)
    logical :: valid,do_transport
    integer :: heat_supplied,heat_min,heat_max,heat_status
    ok=.false.;nv=5+size(species)
    call agree_mpi_thermal_slip_faces(isothermal_slip_faces,comm,thermal_slip_check,ok, &
      faces,transport_flags(2),wall_temperature)
    if(.not.ok) return
    ok=.false.
    mol=.false.;lf=.false.;reconstruction_order=2
    if(present(mol_mode)) mol=mol_mode
    if(present(mol_order)) reconstruction_order=mol_order
    if(present(mol_lax_friedrichs)) lf=mol_lax_friedrichs
    sample_held=.false.
    if(present(auxiliary_samples_internal_energy)) sample_held=auxiliary_samples_internal_energy
    call agree_tree_energy_policy(mol,internal_energy,energy_correction,energy_policy,comm,policy,stored,valid)
    if(.not.valid) return
    heat_supplied=merge(1,0,present(heat_release))
    call MPI_Allreduce(heat_supplied,heat_min,1,MPI_INTEGER,MPI_MIN,comm,heat_status)
    if(heat_status/=MPI_SUCCESS) return
    call MPI_Allreduce(heat_supplied,heat_max,1,MPI_INTEGER,MPI_MAX,comm,heat_status)
    if(heat_status/=MPI_SUCCESS.or.heat_min/=heat_max) return
    call agree_mpi_constant_transport_coefficients(constant_coefficients,comm,ok,species=species)
    if(.not.ok) return
    ok=.false.
    controls=[spacing,dt,force,rtol,atol]
    valid=reconstruction_order>=1.and.reconstruction_order<=2
    valid=valid.and.size(boundary_change)==nv.and.size(source_change)==nv.and.iterations>0
    valid=valid.and.all(ieee_is_finite(controls))
    valid=valid.and.minval(spacing)>0.and.min(dt,rtol,atol)>0.and.valid_reactive_boundaries_3d(faces)
    valid=valid.and..not.(transport_flags(4).and..not.transport_flags(3))
    valid=valid.and.(chemistry.or..not.present(reaction_source))
    valid=valid.and.(chemistry.or..not.present(evaluate_reaction))
    valid=valid.and.(.not.sample_held.or.(stored.and.present(auxiliary_source)))
    valid=valid.and.(present(auxiliary_source).eqv.present(auxiliary_change))
    if(present(auxiliary_change)) valid=valid.and.size(auxiliary_change)==nv
    call agree(valid)
    if(.not.valid) return
    switches=[iterations,merge(1,0,chemistry),merge(1,0,transport_flags),merge(1,0,present(reaction_source)), &
      size(species),size(reactions),merge(1,0,present(evaluate_diffusion)),merge(1,0,present(evaluate_reaction)), &
      merge(1,0,present(auxiliary_source)),merge(1,0,mol),reconstruction_order,merge(1,0,lf), &
      merge(1,0,present(frame)),merge(1,0,present(domain_lower)), &
      merge(1,0,present(auxiliary_samples_internal_energy)),merge(1,0,sample_held)]
    call MPI_Allreduce(switches,low_switches,size(switches),MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(switches,high_switches,size(switches),MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(low_switches/=high_switches)) return
    valid=present(frame).eqv.present(domain_lower)
    call agree(valid)
    if(.not.valid) return
    call MPI_Allreduce(controls,low_controls,9,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high_controls,9,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(low_controls/=high_controls)) return
    if(present(frame)) then
      valid=mol.and.present(domain_lower).and.frame%is_valid()
      if(present(domain_lower)) valid=valid.and.all(ieee_is_finite(domain_lower))
      call agree(valid)
      if(.not.valid) return
      call amr_tree_integer_layout_3d(tree,frame_origin,frame_domain,frame_ratios,valid)
      call agree(valid)
      if(.not.valid) return
      allocate(frame_h(3,size(tree%nodes)))
      frame_h(:,1)=spacing
      do node=2,size(tree%nodes)
        frame_h(:,node)=frame_h(:,tree%nodes(node)%parent)/real(tree%nodes(node)%patch%refinement_ratio,dp)
      enddo
    endif
    depth=-1
    if(present(target_depth)) depth=target_depth
    selection=[depth,merge(1,0,present(coarse_start)),merge(1,0,present(coarse_end)), &
      merge(1,0,present(fractions_start)),merge(1,0,present(fractions_end)), &
      merge(1,0,present(source_start)),merge(1,0,present(source_end)), &
      merge(1,0,present(integrated_flux)),merge(1,0,present(old_predictor)),merge(1,0,present(new_predictor)), &
      merge(1,0,present(internal_energy_integrated_flux)),merge(1,0,present(physical_integrated_flux)), &
      merge(1,0,present(physical_internal_energy_integrated_flux)), &
      merge(1,0,present(internal_energy_pressure_work_change))]
    call MPI_Allreduce(selection,selection_low,size(selection),MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(selection,selection_high,size(selection),MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(selection_low/=selection_high)) return
    if(present(internal_energy_pressure_work_change).and..not.stored) return
    valid=depth>=-1
    if(present(target_depth)) valid=depth>=0.and.present(integrated_flux)
    if(depth>0) valid=valid.and.present(coarse_start)
    valid=valid.and.all(selection(2:7)==selection(2))
    call agree(valid)
    if(.not.valid) return
    if(present(coarse_start)) then
      valid=valid.and.depth>=0
      valid=valid.and.size(fractions_start)==depth.and.size(fractions_end)==depth
    end if
    call agree(valid)
    if(.not.valid) return
    if(present(coarse_start)) then
      valid=all(ieee_is_finite(fractions_start)).and.all(ieee_is_finite(fractions_end))
      valid=valid.and.all(fractions_start>=0).and.all(fractions_end<=1).and.all(fractions_end>=fractions_start)
      call agree(valid)
      if(.not.valid) return
      midpoint=.5_dp*(fractions_start+fractions_end)
    end if
    if((mol.or.stored).and.depth<0) then
      if(present(auxiliary_source)) then
      call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,state,temperature,spacing,dt,faces, &
        iterations,chemistry,transport_flags,force,rtol,atol,boundary_change,source_change,maximum_correction,minimum_theta, &
        comm,ok,reaction_source,inflow_state,inflow_temperature,wall_temperature,wall_velocity,outlet_pressure, &
        inflow_mach,inflow_pressure,difmag,evaluate_diffusion=evaluate_diffusion,evaluate_reaction=evaluate_reaction, &
        constant_coefficients=constant_coefficients,auxiliary_source=ordered_source_adapter,auxiliary_change=auxiliary_change, &
        heat_release=heat_release,mol_mode=mol_mode,mol_order=mol_order,mol_lax_friedrichs=mol_lax_friedrichs, &
        subcycling=.false.,integrated_flux=integrated_flux,old_predictor=old_predictor,new_predictor=new_predictor, &
        internal_energy=internal_energy,energy_correction=energy_correction,energy_policy=energy_policy, &
        internal_energy_integrated_flux=internal_energy_integrated_flux, &
        physical_integrated_flux=physical_integrated_flux, &
        physical_internal_energy_integrated_flux=physical_internal_energy_integrated_flux,frame=frame,domain_lower=domain_lower, &
        internal_energy_pressure_work_change=internal_energy_pressure_work_change, &
        auxiliary_samples_internal_energy=sample_held, &
          isothermal_slip_faces=isothermal_slip_faces)
      else
      call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,state,temperature,spacing,dt,faces, &
        iterations,chemistry,transport_flags,force,rtol,atol,boundary_change,source_change,maximum_correction,minimum_theta, &
        comm,ok,reaction_source,inflow_state,inflow_temperature,wall_temperature,wall_velocity,outlet_pressure, &
        inflow_mach,inflow_pressure,difmag,evaluate_diffusion=evaluate_diffusion,evaluate_reaction=evaluate_reaction, &
        constant_coefficients=constant_coefficients,auxiliary_change=auxiliary_change, &
        heat_release=heat_release,mol_mode=mol_mode,mol_order=mol_order,mol_lax_friedrichs=mol_lax_friedrichs, &
        subcycling=.false.,integrated_flux=integrated_flux,old_predictor=old_predictor,new_predictor=new_predictor, &
        internal_energy=internal_energy,energy_correction=energy_correction,energy_policy=energy_policy, &
        internal_energy_integrated_flux=internal_energy_integrated_flux, &
        physical_integrated_flux=physical_integrated_flux, &
        physical_internal_energy_integrated_flux=physical_internal_energy_integrated_flux,frame=frame,domain_lower=domain_lower, &
        internal_energy_pressure_work_change=internal_energy_pressure_work_change, &
          isothermal_slip_faces=isothermal_slip_faces)
      endif
      return
    endif
    call MPI_Comm_rank(comm,rank,status)
    if(present(coarse_start)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,coarse_start,nv+1+merge(1,0,stored),comm,valid,maximum_depth=depth-1)
      if(.not.valid) return
      call validate_mpi_amr_tree_local_3d(tree,owners,coarse_end,nv+1+merge(1,0,stored),comm,valid,maximum_depth=depth-1)
      if(.not.valid) return
      allocate(transport_start,transport_end)
      allocate(transport_start%fields(size(tree%nodes)),transport_end%fields(size(tree%nodes)))
      valid=.true.
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        if(tree%nodes(node)%level>=depth) then
          if(.not.allocated(coarse_start%fields(node)%values).or..not.allocated(coarse_end%fields(node)%values)) cycle
        endif
        if(any(shape(coarse_start%fields(node)%values)/=[nv+1+merge(1,0,stored),tree%nodes(node)%shape]).or. &
          any(shape(coarse_end%fields(node)%values)/=[nv+1+merge(1,0,stored),tree%nodes(node)%shape])) then
          valid=.false.;exit
        endif
        transport_start%fields(node)%values=coarse_start%fields(node)%values
        transport_end%fields(node)%values=coarse_end%fields(node)%values
      enddo
      call agree(valid)
      if(.not.valid) return
    endif
    do_transport=any(transport_flags(:3));theta=1
    call clone_mpi_amr_tree_inputs_3d(tree,owners,size(species),state,temperature,old,ot,comm,rank,valid)
    if(.not.valid) return
    if(depth>=tree%level_count()) return
    if(depth<0) then
      call average_down_mpi_amr_tree_3d(tree,owners,old,nv,comm,valid)
      if(.not.valid) return
    end if
    energy_delta=0;initial_energy_delta=0;predictor_energy_delta=0;final_energy_delta=0
    pressure_delta=0
    if(present(internal_energy_pressure_work_change).and..not.mol) allocate(pressure_work)
    if(stored) then
      call validate_mpi_amr_tree_local_3d(tree,owners,internal_energy,1,comm,valid)
      if(.not.valid) return
      allocate(old_held,held);old_held=internal_energy
      if(mol.and..not.present(frame).and.policy%update_total_from_internal) then
        call validate_mpi_mol_old_state(species,tree,owners,old,ot,old_held,comm,valid,target_depth=depth)
      else
        call recover(old,ot,valid,old_held,initial_energy_delta)
      endif
      energy_delta=initial_energy_delta;held=old_held
    else
      call recover(old,ot,valid)
    endif
    if(.not.valid) return
    allocate(ir%fields(size(tree%nodes)),stat=ios)
    if(ios==0) then
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        allocate(ir%fields(node)%values(nv-4,n(1),n(2),n(3)),stat=ios)
        if(ios/=0) exit
        ir%fields(node)%values=0
      end do
    end if
    valid=ios==0;call agree(valid)
    if(.not.valid) return
    if(present(heat_release)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,heat_release,1,comm,valid)
      if(.not.valid) return
      allocate(heat);heat=heat_release
    endif
    if(present(heat_release).and..not.chemistry) then
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
        heat%fields(node)%values=0
      enddo
    endif
    if(present(reaction_source)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,reaction_source,nv-4,comm,valid)
      if(.not.valid) return
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        if(.not.all(ieee_is_finite(reaction_source%fields(node)%values))) valid=.false.
      end do
      call agree(valid)
      if(.not.valid) return
      ir=reaction_source
    end if
    call sources(old,ot,old_ext,old_tr,old_flux,valid,old_aux,0,fractions_start)
    if(.not.valid) return
    call predictor(old_ext,old_tr,old_aux,op,valid)
    if(.not.valid) return
    if(mol) then
      if(stored) then
      call compute_mol_mpi_amr_tree_hydro_3d(tree,owners,species,old,ot,spacing,dt,faces,mol_old_hydro,mol_old_flux,comm,valid, &
        inflow_state,inflow_temperature,outlet_pressure,inflow_mach,inflow_pressure, &
        target_depth,coarse_start,coarse_end,fractions_start,reconstruction_order,lf, &
        internal_energy=old_held,internal_energy_rate=old_arate,internal_energy_integrated_flux=old_aflux, &
        frame=frame,domain_lower=domain_lower)
      else
      call compute_mol_mpi_amr_tree_hydro_3d(tree,owners,species,old,ot,spacing,dt,faces,mol_old_hydro,mol_old_flux,comm,valid, &
        inflow_state,inflow_temperature,outlet_pressure,inflow_mach,inflow_pressure, &
        target_depth,coarse_start,coarse_end,fractions_start,reconstruction_order,lf,frame=frame,domain_lower=domain_lower)
      endif
      if(.not.valid) return
    endif
    new_ext=old_ext;new_tr=old_tr;new_aux=old_aux;candidate=old;ct=ot
    do iteration=1,merge(iterations,1,chemistry.or..not.mol)
      if(.not.mol.or.iteration==1) then
      call average_rates()
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
        rates%fields(node)%values(6:,:,:,:)=rates%fields(node)%values(6:,:,:,:)+ir%fields(node)%values(:nv-5,:,:,:)
        rates%fields(node)%values(5,:,:,:)=rates%fields(node)%values(5,:,:,:)+ir%fields(node)%values(nv-4,:,:,:)
      end do
      if(mol) then
        hydro=mol_old_hydro
      else
      if(stored) then
        call compute_sdc_mpi_amr_tree_hydro_3d(tree,owners,species,old,ot,rates,spacing,dt,faces,hydro,flux,comm,valid, &
          inflow_state,inflow_temperature,outlet_pressure,inflow_mach,inflow_pressure,difmag, &
          target_depth,coarse_start,coarse_end,fractions_start,source_start,source_end,midpoint, &
          internal_energy=old_held,internal_energy_rate=old_arate,internal_energy_integrated_flux=old_aflux, &
          internal_energy_pressure_work=pressure_work)
      else
      call compute_sdc_mpi_amr_tree_hydro_3d(tree,owners,species,old,ot,rates,spacing,dt,faces,hydro,flux,comm,valid, &
        inflow_state,inflow_temperature,outlet_pressure,inflow_mach,inflow_pressure,difmag, &
        target_depth,coarse_start,coarse_end,fractions_start,source_start,source_end,midpoint)
      endif
      if(.not.valid) return
      endif
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
        candidate%fields(node)%values=old%fields(node)%values+dt*(rates%fields(node)%values+hydro%fields(node)%values)
      end do
      if(stored) then
        do node=1,size(tree%nodes)
          if(owners(node)/=rank) cycle
          if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
          held%fields(node)%values=old_held%fields(node)%values+dt*old_arate%fields(node)%values
        enddo
        call recover(candidate,ct,valid,held,predictor_energy_delta,predictor_reset)
        predictor_gas=candidate;predictor_held=held
      else
        call recover(candidate,ct,valid)
      endif
      if(.not.valid) return
      endif
      call sources(candidate,ct,new_ext,new_tr,new_flux,valid,new_aux,iteration,fractions_end)
      if(.not.valid) return
      if(mol) then
        if(stored) then
        call compute_mol_mpi_amr_tree_hydro_3d(tree,owners,species,candidate,ct,spacing,dt,faces, &
          mol_new_hydro,mol_new_flux,comm,valid, &
          inflow_state,inflow_temperature,outlet_pressure,inflow_mach,inflow_pressure, &
          target_depth,coarse_start,coarse_end,fractions_end,reconstruction_order,lf, &
          internal_energy=held,internal_energy_rate=new_arate,internal_energy_integrated_flux=new_aflux, &
          frame=frame,domain_lower=domain_lower)
        else
        call compute_mol_mpi_amr_tree_hydro_3d(tree,owners,species,candidate,ct,spacing,dt,faces, &
          mol_new_hydro,mol_new_flux,comm,valid, &
          inflow_state,inflow_temperature,outlet_pressure,inflow_mach,inflow_pressure, &
          target_depth,coarse_start,coarse_end,fractions_end,reconstruction_order,lf,frame=frame,domain_lower=domain_lower)
        endif
        if(.not.valid) return
        hydro=mol_old_hydro;flux=mol_old_flux
        do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
          if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
          hydro%fields(node)%values=.5_dp*(mol_old_hydro%fields(node)%values+mol_new_hydro%fields(node)%values)
          flux(node)%x=.5_dp*(mol_old_flux(node)%x+mol_new_flux(node)%x)
          flux(node)%y=.5_dp*(mol_old_flux(node)%y+mol_new_flux(node)%y)
          flux(node)%z=.5_dp*(mol_old_flux(node)%z+mol_new_flux(node)%z)
          if(iteration>1) new_ext%fields(node)%values=0
        enddo
      endif
      call average_rates()
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
        rates%fields(node)%values=rates%fields(node)%values+hydro%fields(node)%values
      end do
      if(mol.and.chemistry.and.iteration==1) then
        do node=1,size(tree%nodes)
          if(owners(node)/=rank) cycle
          if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
          ! Match the first pinned-PeleC MOL corrector: form the updated state
          ! before recovering its nonreacting source. Later iterations average
          ! source terms directly; these forms differ at floating-point precision.
          rates%fields(node)%values=.5_dp*candidate%fields(node)%values+.5_dp*old%fields(node)%values
          rates%fields(node)%values=rates%fields(node)%values+(.5_dp*dt)*( &
            mol_new_hydro%fields(node)%values+new_ext%fields(node)%values+ &
            new_tr%fields(node)%values+new_aux%fields(node)%values)
          rates%fields(node)%values(6:,:,:,:)=rates%fields(node)%values(6:,:,:,:)+ &
            (.5_dp*dt)*ir%fields(node)%values(:nv-5,:,:,:)
          rates%fields(node)%values(5,:,:,:)=rates%fields(node)%values(5,:,:,:)+ &
            (.5_dp*dt)*ir%fields(node)%values(nv-4,:,:,:)
          rates%fields(node)%values=(1._dp/dt)*rates%fields(node)%values+(-1._dp/dt)*old%fields(node)%values
          rates%fields(node)%values(6:,:,:,:)=rates%fields(node)%values(6:,:,:,:)-ir%fields(node)%values(:nv-5,:,:,:)
          rates%fields(node)%values(5,:,:,:)=rates%fields(node)%values(5,:,:,:)-ir%fields(node)%values(nv-4,:,:,:)
        enddo
      endif
      if(stored.and.mol.and.iteration==1.and..not.chemistry) then
        do node=1,size(tree%nodes)
          if(owners(node)/=rank) cycle
          if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
          rates%fields(node)%values(5,:,:,:)=rates%fields(node)%values(5,:,:,:)+ &
            .5_dp*predictor_reset%fields(node)%values(1,:,:,:)/dt
        enddo
      endif
      if(chemistry) then
        call react_sdc_mpi_amr_tree_3d(tree,owners,species,reactions,old,ot,rates,dt,rtol,atol,candidate,ct,ir,comm,valid, &
          target_depth,evaluate_reaction,heat_release=heat,internal_energy=held, &
          frame=frame,domain_lower=domain_lower,root_spacing=spacing,old_internal_energy=old_held)
        if(.not.valid) return
        if(stored) call recover(candidate,ct,valid,held,final_energy_delta)
        if(.not.valid) return
      else
        do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
          if(stored.and.mol) then
            candidate%fields(node)%values=.5_dp*(old%fields(node)%values+predictor_gas%fields(node)%values)+ &
              .5_dp*dt*(mol_new_hydro%fields(node)%values+new_ext%fields(node)%values+ &
              new_tr%fields(node)%values+new_aux%fields(node)%values)
            held%fields(node)%values=.5_dp*(old_held%fields(node)%values+predictor_held%fields(node)%values)+ &
              .5_dp*dt*new_arate%fields(node)%values
          else
          candidate%fields(node)%values=old%fields(node)%values+dt*rates%fields(node)%values
          if(stored) held%fields(node)%values=old_held%fields(node)%values+dt*old_arate%fields(node)%values
          endif
        end do
        if(stored) then
          call recover(candidate,ct,valid,held,final_energy_delta)
        else
          call recover(candidate,ct,valid)
        endif
        if(.not.valid) return
      end if
    end do
    if(allocated(pressure_work)) then
      call inventory(pressure_work,pressure_integral,valid)
      if(.not.valid) return
      if(dt>1) then
        if(abs(pressure_integral(1))>huge(1._dp)/dt) return
      endif
      pressure_delta=-dt*pressure_integral(1)
      if(.not.ieee_is_finite(pressure_delta)) return
    endif
    if(stored) then
      energy_delta=energy_delta+final_energy_delta
      if(mol.and.(.not.chemistry.or.iterations==1)) energy_delta=energy_delta+.5_dp*predictor_energy_delta
      aflux=old_aflux;physical_aflux=old_aflux
      if(mol) then
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
        physical_aflux(node)%x=.5_dp*(old_aflux(node)%x+new_aflux(node)%x)
        physical_aflux(node)%y=.5_dp*(old_aflux(node)%y+new_aflux(node)%y)
        physical_aflux(node)%z=.5_dp*(old_aflux(node)%z+new_aflux(node)%z)
        if(.not.chemistry.and.iterations>1) then
          aflux(node)%x=.5_dp*old_aflux(node)%x
          aflux(node)%y=.5_dp*old_aflux(node)%y
          aflux(node)%z=.5_dp*old_aflux(node)%z
        else
          aflux(node)%x=.5_dp*(old_aflux(node)%x+new_aflux(node)%x)
          aflux(node)%y=.5_dp*(old_aflux(node)%y+new_aflux(node)%y)
          aflux(node)%z=.5_dp*(old_aflux(node)%z+new_aflux(node)%z)
        endif
      enddo
      endif
    endif
    if(do_transport.or.present(evaluate_diffusion)) then
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
        flux(node)%x=flux(node)%x+.5_dp*(old_flux(node)%x+new_flux(node)%x)
        flux(node)%y=flux(node)%y+.5_dp*(old_flux(node)%y+new_flux(node)%y)
        flux(node)%z=flux(node)%z+.5_dp*(old_flux(node)%z+new_flux(node)%z)
      end do
    end if
    registration=flux
    if(mol.and..not.chemistry.and.iterations>1) then
      ! The fixed reference registers only the old half when no reacting
      ! iteration exists to register the final new half. Keep physical flux
      ! for boundary accounting and the actual register for synchronization.
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
        registration(node)%x=.5_dp*mol_old_flux(node)%x
        registration(node)%y=.5_dp*mol_old_flux(node)%y
        registration(node)%z=.5_dp*mol_old_flux(node)%z
        if(do_transport.or.present(evaluate_diffusion)) then
          registration(node)%x=registration(node)%x+.5_dp*old_flux(node)%x
          registration(node)%y=registration(node)%y+.5_dp*old_flux(node)%y
          registration(node)%z=registration(node)%z+.5_dp*old_flux(node)%z
        endif
      enddo
    endif
    budget=0;largest=0
    if(rank==owners(1).and.depth<=0) then
    n=tree%nodes(1)%shape
    ! Subtract opposing faces before summing to retain small net impulses.
    do c=1,nv
      budget(c)=sum(flux(1)%x(c,1,:,:)-flux(1)%x(c,n(1)+1,:,:))*spacing(2)*spacing(3)+ &
        sum(flux(1)%y(c,:,1,:)-flux(1)%y(c,:,n(2)+1,:))*spacing(1)*spacing(3)+ &
        sum(flux(1)%z(c,:,:,1)-flux(1)%z(c,:,:,n(3)+1))*spacing(1)*spacing(2)
    end do
    end if
    call MPI_Bcast(budget,nv,MPI_DOUBLE_PRECISION,owners(1),comm,status)
    if(status/=MPI_SUCCESS) return
    if(depth<0) then
    if(stored) then
    call reflux_mpi_amr_tree_3d(tree,owners,candidate,nv,registration,spacing,faces,correction,largest,comm,valid, &
      first_species=6,last_species=nv,synchronized_flux=effective_mol_flux,auxiliary_state=held,auxiliary_flux=aflux, &
      auxiliary_synchronized_flux=effective_mol_aux)
    else if(mol) then
    call reflux_mpi_amr_tree_3d(tree,owners,candidate,nv,registration,spacing,faces,correction,largest,comm,valid, &
      first_species=6,last_species=nv,synchronized_flux=effective_mol_flux)
    else
    call reflux_mpi_amr_tree_3d(tree,owners,candidate,nv,registration,spacing,faces,correction,largest,comm,valid)
    endif
    if(.not.valid) return
    if(mol.or.stored) then
      ! Add the synchronization change to the physical two-half-step flux,
      ! retaining the pinned distinction between that flux and its register.
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        flux(node)%x=flux(node)%x+(effective_mol_flux(node)%x-registration(node)%x)
        flux(node)%y=flux(node)%y+(effective_mol_flux(node)%y-registration(node)%y)
        flux(node)%z=flux(node)%z+(effective_mol_flux(node)%z-registration(node)%z)
        if(stored) then
          physical_aflux(node)%x=physical_aflux(node)%x+(effective_mol_aux(node)%x-aflux(node)%x)
          physical_aflux(node)%y=physical_aflux(node)%y+(effective_mol_aux(node)%y-aflux(node)%y)
          physical_aflux(node)%z=physical_aflux(node)%z+(effective_mol_aux(node)%z-aflux(node)%z)
        endif
      enddo
      registration=effective_mol_flux
      if(stored) aflux=effective_mol_aux
    endif
    budget=budget+correction
    call recover(candidate,ct,valid)
    if(.not.valid) return
    end if
    call inventory(old_ext,first,valid)
    if(.not.valid) return
    call inventory(new_ext,last,valid)
    if(.not.valid) return
    volume=.5_dp*dt*(first+last)
    if(.not.all(ieee_is_finite(budget)).or..not.all(ieee_is_finite(volume))) return
    if(present(auxiliary_source)) then
      call inventory(old_aux,first,valid)
      if(.not.valid) return
      call inventory(new_aux,last,valid)
      if(.not.valid) return
      auxiliary=.5_dp*dt*(first+last)
      if(.not.all(ieee_is_finite(auxiliary))) return
    endif
    if(present(new_predictor)) then
      call predictor(new_ext,new_tr,new_aux,np,valid)
      if(.not.valid) return
    end if
    if(present(reaction_source).and.depth<0) then
      call average_down_mpi_amr_tree_3d(tree,owners,ir,nv-4,comm,valid)
      if(.not.valid) return
    endif
    if(present(heat_release).and.depth<0) then
      call average_down_mpi_amr_tree_3d(tree,owners,heat,1,comm,valid)
      if(.not.valid) return
    endif
    if(mol) then
      valid=.true.
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
        if(any(candidate%fields(node)%values(6:,:,:,:)<0._dp)) valid=.false.
      enddo
      call agree(valid)
      if(.not.valid) return
    endif
    if(stored) then
      valid=ieee_is_finite(energy_correction+energy_delta);call agree(valid)
      if(.not.valid) return
      internal_energy=held;energy_correction=energy_correction+energy_delta
      if(present(internal_energy_integrated_flux)) internal_energy_integrated_flux=aflux
      if(present(physical_internal_energy_integrated_flux)) &
        physical_internal_energy_integrated_flux=physical_aflux
    endif
    if(present(internal_energy_pressure_work_change)) internal_energy_pressure_work_change=pressure_delta
    state=candidate;temperature=ct;boundary_change=budget;source_change=volume
    maximum_correction=largest;minimum_theta=theta
    if(present(auxiliary_change)) auxiliary_change=auxiliary
    if(present(heat_release)) heat_release=heat
    if(present(reaction_source)) reaction_source=ir
    if(present(integrated_flux)) integrated_flux=registration
    if(present(physical_integrated_flux)) physical_integrated_flux=flux
    if(present(old_predictor)) old_predictor=op
    if(present(new_predictor)) new_predictor=np
    ok=.true.
  contains
    subroutine ordered_source_adapter(stage_state,stage_temperature,interval_dt,stage,depth,time,rate,source_comm,success)
      type(amr_tree_state_3d), intent(in) :: stage_state,stage_temperature
      type(amr_tree_state_3d), intent(inout) :: rate
      real(dp), intent(in) :: interval_dt,time
      integer, intent(in) :: stage,depth
      type(MPI_Comm), intent(in) :: source_comm
      logical, intent(out) :: success
      call auxiliary_source(stage_state,stage_temperature,interval_dt,stage,depth,rate,source_comm,success)
    end subroutine ordered_source_adapter
    subroutine predictor(external,diffusive,aux,result,good)
      type(amr_tree_state_3d), intent(in) :: external,diffusive,aux
      type(amr_tree_state_3d), intent(out) :: result
      logical, intent(out) :: good
      integer :: peer
      result=external;good=.true.
      do peer=1,size(tree%nodes)
        if(owners(peer)/=rank) cycle
        if(depth>=0.and.tree%nodes(peer)%level/=depth) cycle
        result%fields(peer)%values=result%fields(peer)%values+diffusive%fields(peer)%values
        if(present(auxiliary_source)) result%fields(peer)%values=result%fields(peer)%values+aux%fields(peer)%values
        result%fields(peer)%values(6:,:,:,:)=result%fields(peer)%values(6:,:,:,:)+ir%fields(peer)%values(:nv-5,:,:,:)
        result%fields(peer)%values(5,:,:,:)=result%fields(peer)%values(5,:,:,:)+ir%fields(peer)%values(nv-4,:,:,:)
        if(.not.all(ieee_is_finite(result%fields(peer)%values))) good=.false.
      end do
      call agree(good)
    end subroutine predictor
    subroutine average_rates()
      integer :: peer
      rates=old_ext
      do peer=1,size(tree%nodes)
        if(owners(peer)/=rank) cycle
        rates%fields(peer)%values=.5_dp*new_ext%fields(peer)%values+.5_dp*old_ext%fields(peer)%values
        rates%fields(peer)%values=rates%fields(peer)%values+.5_dp*new_tr%fields(peer)%values
        rates%fields(peer)%values=rates%fields(peer)%values+.5_dp*old_tr%fields(peer)%values
        if(present(auxiliary_source)) rates%fields(peer)%values=rates%fields(peer)%values+ &
          .5_dp*(old_aux%fields(peer)%values+new_aux%fields(peer)%values)
      end do
    end subroutine
    subroutine recover(u,t,good,held_energy,change,delta)
      type(amr_tree_state_3d), intent(inout) :: u,t
      type(amr_tree_state_3d), optional, intent(inout) :: held_energy
      real(dp), optional, intent(out) :: change
      type(amr_tree_state_3d), optional, intent(out) :: delta
      logical, intent(out) :: good
      real(dp), allocatable :: recovered(:,:,:),potential(:,:,:)
      real(dp) :: position(3)
      integer :: peer,dims(3),i,j,k
      if(present(held_energy)) then
        call reset_mpi_tree_energy(species,tree,owners,u,t,held_energy,spacing,policy,change,comm,good, &
          target_depth=depth,delta=delta,frame=frame,domain_lower=domain_lower)
      else if(present(frame)) then
        good=.true.
        do peer=1,size(tree%nodes)
          if(owners(peer)/=rank) cycle
          if(depth>=0.and.tree%nodes(peer)%level/=depth) cycle
          dims=tree%nodes(peer)%shape
          allocate(recovered(dims(1),dims(2),dims(3)),potential(dims(1),dims(2),dims(3)))
          do k=1,dims(3)
          do j=1,dims(2)
          do i=1,dims(1)
            position=domain_lower+(real(frame_origin(:,peer),dp)+real([i,j,k],dp)-.5_dp)*frame_h(:,peer)
            call rotating_potential_3d(frame,position,potential(i,j,k),good)
            if(.not.good) exit
          enddo
          if(.not.good) exit
          enddo
          if(.not.good) exit
          enddo
          if(good) then
            call recover_reactive_temperatures_3d(species,u%fields(peer)%values,t%fields(peer)%values(1,:,:,:), &
              dims(1),dims(2),dims(3),recovered,good,rotating_potential=potential)
          endif
          if(.not.good) exit
          t%fields(peer)%values(1,:,:,:)=recovered
          deallocate(recovered,potential)
        enddo
        call agree(good)
      else
        call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,u,t,comm,good,target_depth)
      endif
    end subroutine
    subroutine sources(u,t,external,diffusive,face_flux,good,aux,stage,fractions)
      type(amr_tree_state_3d), intent(in) :: u,t
      type(amr_tree_state_3d), intent(inout) :: external,diffusive,aux
      type(amr_tree_flux_3d), allocatable, intent(inout) :: face_flux(:)
      logical, intent(out) :: good
      real(dp), optional, intent(in) :: fractions(:)
      integer :: peer,cells(3),stage,i,j,k
      real(dp) :: position(3),rotation_source(3)
      type(amr_tree_state_3d), allocatable :: stage_held
      type(amr_tree_state_3d) :: extra,sampling
      integer :: allocation_status
      type(amr_tree_flux_3d), allocatable :: extra_flux(:)
      external=u;diffusive=u;good=.true.
      do peer=1,size(tree%nodes)
        if(owners(peer)/=rank) cycle
        external%fields(peer)%values=0;diffusive%fields(peer)%values=0
        if(depth>=0.and.tree%nodes(peer)%level/=depth) cycle
        call compute_reactive_external_source_3d(u%fields(peer)%values,force,external%fields(peer)%values,good)
        if(.not.good) exit
        if(present(frame)) then
          cells=tree%nodes(peer)%shape
          do k=1,cells(3)
          do j=1,cells(2)
          do i=1,cells(1)
            position=domain_lower+(real(frame_origin(:,peer),dp)+real([i,j,k],dp)-.5_dp)*frame_h(:,peer)
            call rotating_momentum_source_3d(frame,position,u%fields(peer)%values(1,i,j,k), &
              u%fields(peer)%values(2:4,i,j,k),rotation_source,good)
            if(.not.good) exit
            external%fields(peer)%values(2:4,i,j,k)=external%fields(peer)%values(2:4,i,j,k)+rotation_source
          enddo
          if(.not.good) exit
          enddo
          if(.not.good) exit
          enddo
          if(.not.good) exit
        endif
        diffusive%fields(peer)%values=0
      end do
      call agree(good)
      if(.not.good) return
      aux=diffusive
      if(present(auxiliary_source)) then
        ! Callback rates use the same conserved energy convention as u.
        ! Pinned reacting MOL reevaluates the spatial operator after stage 1,
        ! without rebuilding the new external particle/film source family.
        if(.not.(mol.and.stage>1)) then
          if(sample_held) then
            sampling=u;good=.true.
            do peer=1,size(tree%nodes)
              if(owners(peer)/=rank) cycle
              cells=tree%nodes(peer)%shape
              deallocate(sampling%fields(peer)%values)
              allocate(sampling%fields(peer)%values(nv+1,cells(1),cells(2),cells(3)),stat=allocation_status)
              if(allocation_status/=0) then
                good=.false.;exit
              endif
              sampling%fields(peer)%values(:nv,:,:,:)=u%fields(peer)%values
              if(stage==0) then
                sampling%fields(peer)%values(nv+1:nv+1,:,:,:)=old_held%fields(peer)%values
              else
                sampling%fields(peer)%values(nv+1:nv+1,:,:,:)=held%fields(peer)%values
              endif
              good=all(ieee_is_finite(sampling%fields(peer)%values))
              if(.not.good) exit
            enddo
            call agree(good)
            if(.not.good) return
            call auxiliary_source(sampling,t,dt,stage,depth,aux,comm,good)
          else
          call auxiliary_source(u,t,dt,stage,depth,aux,comm,good)
          endif
          call agree(good)
          if(.not.good) return
        endif
        call validate_mpi_amr_tree_local_3d(tree,owners,aux,nv,comm,good)
        if(.not.good) return
        do peer=1,size(tree%nodes)
          if(owners(peer)/=rank) cycle
          if(depth>=0.and.tree%nodes(peer)%level/=depth) aux%fields(peer)%values=0
          if(.not.all(ieee_is_finite(aux%fields(peer)%values))) good=.false.
        enddo
        call agree(good)
        if(.not.good) return
      endif
      if(stored) then
        allocate(stage_held)
        if(stage==0) then
          stage_held=old_held
        else
          stage_held=held
        endif
      endif
      if(do_transport) then
        stage_theta=1
        call compute_sdc_mpi_amr_tree_transport_3d(tree,owners,species,transport,u,t,spacing,dt,faces,transport_flags, &
          diffusive,face_flux,stage_theta,comm,good,inflow_state,inflow_temperature,wall_temperature,wall_velocity, &
          outlet_pressure,inflow_mach,inflow_pressure,target_depth,transport_start,transport_end,fractions, &
          constant_coefficients=constant_coefficients,frame=frame,domain_lower=domain_lower, &
          internal_energy=stage_held,total_energy_primitives=mol.and..not.present(frame), &
          isothermal_slip_faces=isothermal_slip_faces)
        if(.not.good) return
        theta=min(theta,stage_theta)
      end if
      if(present(evaluate_diffusion)) then
        extra=u
        call evaluate_diffusion(u,t,depth,fractions,transport_start,transport_end,extra,extra_flux,good,internal_energy=stage_held)
        call agree(good)
        if(.not.good) return
        call validate_mpi_amr_tree_local_3d(tree,owners,extra,nv,comm,good)
        if(.not.good) return
        good=allocated(extra_flux)
        if(good) good=size(extra_flux)==size(tree%nodes)
        call agree(good)
        if(.not.good) return
        do peer=1,size(tree%nodes)
          if(owners(peer)/=rank) cycle
          if(depth>=0.and.tree%nodes(peer)%level/=depth) cycle
          cells=tree%nodes(peer)%shape
          good=allocated(extra_flux(peer)%x).and.allocated(extra_flux(peer)%y).and.allocated(extra_flux(peer)%z)
          if(.not.good) exit
          good=all(shape(extra_flux(peer)%x)==[nv,cells(1)+1,cells(2),cells(3)]).and. &
            all(shape(extra_flux(peer)%y)==[nv,cells(1),cells(2)+1,cells(3)]).and. &
            all(shape(extra_flux(peer)%z)==[nv,cells(1),cells(2),cells(3)+1])
          if(.not.good) exit
          good=all(ieee_is_finite(extra%fields(peer)%values)).and. &
            all(ieee_is_finite(extra_flux(peer)%x)).and.all(ieee_is_finite(extra_flux(peer)%y)).and. &
            all(ieee_is_finite(extra_flux(peer)%z))
          if(.not.good) exit
          if(present(frame)) then
            good=all(extra%fields(peer)%values(1,:,:,:)==0).and.all(extra_flux(peer)%x(1,:,:,:)==0).and. &
              all(extra_flux(peer)%y(1,:,:,:)==0).and.all(extra_flux(peer)%z(1,:,:,:)==0)
            if(.not.good) exit
          endif
          diffusive%fields(peer)%values=diffusive%fields(peer)%values+extra%fields(peer)%values
          extra_flux(peer)%x=dt*extra_flux(peer)%x;extra_flux(peer)%y=dt*extra_flux(peer)%y
          extra_flux(peer)%z=dt*extra_flux(peer)%z
          if(do_transport) then
            face_flux(peer)%x=face_flux(peer)%x+extra_flux(peer)%x
            face_flux(peer)%y=face_flux(peer)%y+extra_flux(peer)%y
            face_flux(peer)%z=face_flux(peer)%z+extra_flux(peer)%z
          end if
        end do
        call agree(good)
        if(.not.good) return
        if(.not.do_transport) face_flux=extra_flux
      end if
      good=.true.
    end subroutine
    subroutine agree(condition)
      logical, intent(inout) :: condition
      logical :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      condition=.false.
      if(ierr==MPI_SUCCESS) condition=result
    end subroutine
    subroutine inventory(fields,result,good)
      type(amr_tree_state_3d), intent(in) :: fields
      real(dp), intent(out) :: result(:)
      logical, intent(out) :: good
      real(dp) :: accum(size(result),2),term(size(result)),updated(size(result)),vol(size(tree%nodes)),ratio
      integer :: peer,parent,ii,jj,kk,dims(3),ierr
      good=.false.;accum=0;vol(1)=product(spacing)
      if(.not.ieee_is_finite(vol(1)).or.vol(1)<=0) return
      do peer=1,size(tree%nodes)
        if(peer>1) then
          parent=tree%nodes(peer)%parent;ratio=real(tree%nodes(peer)%patch%refinement_ratio,dp)
          vol(peer)=vol(parent)/ratio/ratio/ratio
          if(vol(peer)<=0) return
        end if
        if(owners(peer)==rank) then
          dims=tree%nodes(peer)%shape
          do kk=1,dims(3)
          do jj=1,dims(2)
          do ii=1,dims(1)
            if(covered(tree,peer,[ii,jj,kk])) cycle
            term=vol(peer)*fields%fields(peer)%values(:,ii,jj,kk)-accum(:,2)
            updated=accum(:,1)+term
            accum(:,2)=(updated-accum(:,1))-term;accum(:,1)=updated
          end do
          end do
          end do
        end if
        call MPI_Bcast(accum,size(accum),MPI_DOUBLE_PRECISION,owners(peer),comm,ierr)
        if(ierr/=MPI_SUCCESS) return
      end do
      if(.not.all(ieee_is_finite(accum))) return
      result=accum(:,1);good=.true.
    end subroutine
  end subroutine advance_mpi_amr_tree_sdc_3d

  subroutine validate_mpi_mol_old_state(species,tree,owners,state,temperature,held,comm,ok,target_depth)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(amr_tree_state_3d), intent(in) :: state,temperature,held
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer, optional, intent(in) :: target_depth
    real(dp), allocatable :: recovered(:,:,:)
    integer :: peer,n(3),rank,status
    logical :: local_ok,global_ok
    logical :: valid
    ok=.false.
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    call validate_local()
    call MPI_Allreduce(local_ok,global_ok,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status==MPI_SUCCESS) ok=global_ok
  contains
    subroutine validate_local()
    local_ok=.false.
    do peer=1,size(tree%nodes)
      if(owners(peer)/=rank) cycle
      if(present(target_depth)) then
        if(target_depth>=0.and.tree%nodes(peer)%level/=target_depth) cycle
      endif
      n=tree%nodes(peer)%shape
      if(.not.all(ieee_is_finite(state%fields(peer)%values))) return
      if(.not.all(ieee_is_finite(temperature%fields(peer)%values))) return
      if(any(temperature%fields(peer)%values<=0)) return
      if(any(state%fields(peer)%values(6:,:,:,:)<0)) return
      allocate(recovered(n(1),n(2),n(3)))
      call recover_reactive_temperatures_3d(species,state%fields(peer)%values,temperature%fields(peer)%values(1,:,:,:), &
        n(1),n(2),n(3),recovered,valid)
      if(.not.valid) return
      call recover_reactive_temperatures_3d(species,state%fields(peer)%values,temperature%fields(peer)%values(1,:,:,:), &
        n(1),n(2),n(3),recovered,valid,internal_energy_density=held%fields(peer)%values(1,:,:,:))
      if(.not.valid) return
      deallocate(recovered)
    enddo
    local_ok=.true.
    end subroutine validate_local
  end subroutine validate_mpi_mol_old_state

  subroutine agree_tree_energy_policy(mol,internal_energy,energy_correction,energy_policy,comm,policy,stored,ok)
    logical, intent(in) :: mol
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy
    real(dp), optional, intent(in) :: energy_correction
    type(dual_energy_policy), optional, intent(in) :: energy_policy
    type(MPI_Comm), intent(in) :: comm
    type(dual_energy_policy), intent(out) :: policy
    logical, intent(out) :: stored,ok
    integer :: flags(7),lo(7),hi(7),status
    real(dp) :: values(2),low(2),high(2)
    logical :: valid
    policy=dual_energy_policy();stored=present(internal_energy);ok=.false.
    if(present(energy_policy)) policy=energy_policy
    flags=[merge(1,0,stored),merge(1,0,present(energy_correction)),merge(1,0,present(energy_policy)), &
      merge(1,0,mol),merge(1,0,policy%allow_small_energy),merge(1,0,policy%allow_negative_energy), &
      merge(1,0,policy%update_total_from_internal)]
    call MPI_Allreduce(flags,lo,size(flags),MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,hi,size(flags),MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(lo/=hi)) return
    if(flags(1)/=flags(2).or.(flags(3)>0.and..not.stored)) return
    ! Both native MOL and PPM/CTU carry an independent stored UEINT family.
    values=[policy%eta2,0._dp]
    if(present(energy_correction)) values(2)=energy_correction
    valid=all(ieee_is_finite(values)).and.policy%eta2>=0
    call MPI_Allreduce(valid,ok,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS) ok=.false.
    if(.not.ok) return
    ok=.false.
    call MPI_Allreduce(values,low,size(values),MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(values,high,size(values),MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    ok=.true.
  end subroutine agree_tree_energy_policy

  ! Close a caller's physical gas correction with the same native dual-energy
  ! policy. Covered averages and the numerical prefix publish atomically.
  subroutine close_mpi_amr_tree_held_energy_3d(species,tree,owners,state,temperature,internal_energy,spacing, &
      policy,prefix,comm,ok,frame,domain_lower)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature,internal_energy
    real(dp), intent(in) :: spacing(3)
    type(dual_energy_policy), intent(in) :: policy
    real(dp), intent(inout) :: prefix
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3)
    type(amr_tree_state_3d) :: gas,temps,held
    type(dual_energy_policy) :: agreed_policy
    real(dp) :: change,next_prefix,values(10),low(10),high(10)
    integer :: flags(2),lo(2),hi(2),status,nv
    logical :: valid,ready,stored
    ok=.false.;nv=5+size(species)
    call agree_tree_energy_policy(.true.,internal_energy,prefix,policy,comm,agreed_policy,stored,ready)
    if(.not.ready) return
    flags=[merge(1,0,present(frame)),merge(1,0,present(domain_lower))]
    call MPI_Allreduce(flags,lo,2,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,hi,2,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(lo/=hi).or.flags(1)/=flags(2)) return
    values=0;values(:3)=spacing
    if(present(frame)) then
      values(4:6)=domain_lower;values(7)=frame%angular_speed;values(8:10)=frame%axis_point
      flags(1)=frame%axis
    endif
    valid=all(ieee_is_finite(values)).and.all(spacing>0)
    if(present(frame)) valid=valid.and.frame%is_valid()
    call MPI_Allreduce(valid,ready,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.ready) return
    call MPI_Allreduce(flags,lo,2,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,hi,2,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(lo/=hi)) return
    call MPI_Allreduce(values,low,10,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(values,high,10,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nv,comm,ready)
    if(.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,ready)
    if(.not.ready) return
    call validate_mpi_amr_tree_local_3d(tree,owners,internal_energy,1,comm,ready)
    if(.not.ready) return
    gas=state;temps=temperature;held=internal_energy
    call reset_mpi_tree_energy(species,tree,owners,gas,temps,held,spacing,agreed_policy,change,comm,ready, &
      frame=frame,domain_lower=domain_lower)
    if(.not.ready) return
    valid=ieee_is_finite(change)
    if(valid.and.change>0) valid=prefix<=huge(1._dp)-change
    if(valid.and.change<0) valid=prefix>=-huge(1._dp)-change
    call MPI_Allreduce(valid,ready,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS.or..not.ready) return
    next_prefix=prefix+change
    call average_down_mpi_amr_tree_3d(tree,owners,gas,nv,comm,ready)
    if(.not.ready) return
    call average_down_mpi_amr_tree_3d(tree,owners,held,1,comm,ready)
    if(.not.ready) return
    call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,gas,temps,comm,ready,internal_energy=held)
    if(.not.ready) return
    call move_alloc(gas%fields,state%fields);call move_alloc(temps%fields,temperature%fields)
    call move_alloc(held%fields,internal_energy%fields);prefix=next_prefix;ok=.true.
  end subroutine close_mpi_amr_tree_held_energy_3d

  subroutine reset_mpi_tree_energy(species,tree,owners,state,temperature,internal_energy,spacing,policy,change,comm,ok, &
      target_depth,minimum_depth,delta,frame,domain_lower)
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature,internal_energy
    real(dp), intent(in) :: spacing(3)
    type(dual_energy_policy), intent(in) :: policy
    real(dp), optional, intent(out) :: change
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer, optional, intent(in) :: target_depth,minimum_depth
    type(amr_tree_state_3d), optional, intent(out) :: delta
    type(amr_tree_state_3d) :: corrections
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: recovered(:,:,:)
    real(dp) :: cell_change,integral(1),position(3),potential,h(3)
    integer :: peer,n(3),i,j,k,depth,floor_depth,rank,status,ios,ancestor
    logical :: valid
    ok=.false.;depth=-1;floor_depth=0
    if(present(target_depth)) depth=target_depth
    if(present(minimum_depth)) floor_depth=minimum_depth
    if(present(frame)) then
      if(.not.present(domain_lower)) return
      call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
      if(.not.valid) return
    endif
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    corrections=internal_energy;valid=.true.
    nodes: do peer=1,size(tree%nodes)
      if(owners(peer)/=rank) cycle
      if(.not.all(ieee_is_finite(internal_energy%fields(peer)%values))) then
        valid=.false.;exit nodes
      endif
      corrections%fields(peer)%values=0
      if(depth>=0.and.tree%nodes(peer)%level/=depth) cycle
      if(tree%nodes(peer)%level<floor_depth) cycle
      n=tree%nodes(peer)%shape
      if(present(frame)) then
        h=spacing;ancestor=peer
        do while(ancestor>1)
          h=h/real(tree%nodes(ancestor)%patch%refinement_ratio,dp)
          ancestor=tree%nodes(ancestor)%parent
        enddo
      endif
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        cell_change=0
        if(present(frame)) then
          position=domain_lower+(real(origin(:,peer),dp)+real([i,j,k],dp)-.5_dp)*h
          call rotating_potential_3d(frame,position,potential,valid)
          if(.not.valid) exit nodes
          call reset_nasa7_dual_energy_cell(species,state%fields(peer)%values(:,i,j,k), &
            internal_energy%fields(peer)%values(1,i,j,k),temperature%fields(peer)%values(1,i,j,k), &
            policy,cell_change,valid,potential=potential)
        else
          call reset_nasa7_dual_energy_cell(species,state%fields(peer)%values(:,i,j,k), &
            internal_energy%fields(peer)%values(1,i,j,k),temperature%fields(peer)%values(1,i,j,k), &
            policy,cell_change,valid)
        endif
        if(.not.valid) exit nodes
        corrections%fields(peer)%values(1,i,j,k)=cell_change
      enddo
      enddo
      enddo
      allocate(recovered(n(1),n(2),n(3)),stat=ios)
      if(ios/=0) then
        valid=.false.;exit nodes
      endif
      call recover_reactive_temperatures_3d(species,state%fields(peer)%values,temperature%fields(peer)%values(1,:,:,:), &
        n(1),n(2),n(3),recovered,valid,internal_energy_density=internal_energy%fields(peer)%values(1,:,:,:))
      if(.not.valid) exit nodes
      temperature%fields(peer)%values(1,:,:,:)=recovered
      deallocate(recovered)
    enddo nodes
    call MPI_Allreduce(valid,ok,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS) ok=.false.
    if(.not.ok) return
    call tree_energy_inventory(tree,owners,corrections,spacing,comm,rank,integral,ok)
    if(.not.ok) return
    if(present(change)) change=integral(1)
    if(present(delta)) delta=corrections
  end subroutine reset_mpi_tree_energy

  subroutine tree_energy_inventory(tree,owners,fields,spacing,comm,rank,result,good)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),rank
    real(dp), intent(in) :: spacing(3)
    type(MPI_Comm), intent(in) :: comm
      type(amr_tree_state_3d), intent(in) :: fields
      real(dp), intent(out) :: result(:)
      logical, intent(out) :: good
      real(dp) :: accum(size(result),2),term(size(result)),updated(size(result)),vol(size(tree%nodes)),ratio
      integer :: peer,parent,ii,jj,kk,dims(3),ierr
      good=.false.;accum=0;vol(1)=product(spacing)
      if(.not.ieee_is_finite(vol(1)).or.vol(1)<=0) return
      do peer=1,size(tree%nodes)
        if(peer>1) then
          parent=tree%nodes(peer)%parent;ratio=real(tree%nodes(peer)%patch%refinement_ratio,dp)
          vol(peer)=vol(parent)/ratio/ratio/ratio
          if(vol(peer)<=0) return
        end if
        if(owners(peer)==rank) then
          dims=tree%nodes(peer)%shape
          do kk=1,dims(3)
          do jj=1,dims(2)
          do ii=1,dims(1)
            if(covered(tree,peer,[ii,jj,kk])) cycle
            term=vol(peer)*fields%fields(peer)%values(:,ii,jj,kk)-accum(:,2)
            updated=accum(:,1)+term
            accum(:,2)=(updated-accum(:,1))-term;accum(:,1)=updated
          end do
          end do
          end do
        end if
        call MPI_Bcast(accum,size(accum),MPI_DOUBLE_PRECISION,owners(peer),comm,ierr)
        if(ierr/=MPI_SUCCESS) return
      end do
      if(.not.all(ieee_is_finite(accum))) return
      result=accum(:,1);good=.true.
  end subroutine tree_energy_inventory

  ! R/2 - T/2 - H - T/2 - R/2, over one coarse interval, optionally with recursive flow/transport steps.
  ! Disabled operators are skipped. All stage outputs remain private until
  ! the complete composition succeeds, including the final chemistry solve.
  subroutine advance_mpi_amr_tree_split_3d(tree,owners,species,reactions,transport,state,temperature,spacing,dt,faces, &
      reconstruction,limiter,solver,chemistry,viscosity,conduction,diffusion,barodiffusion,rtol,atol, &
      boundary_change,maximum_correction,minimum_theta,comm,ok,chemistry_integrator, &
      inflow_state,inflow_temperature,wall_temperature,wall_velocity,subcycling,outlet_pressure,inflow_mach,inflow_pressure, &
      difmag, &
      constant_coefficients)
    real(dp), optional, intent(in) :: constant_coefficients(:)
    logical :: constant_agreed
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: spacing(3),dt,rtol,atol
    character(len=*), intent(in) :: faces(6),reconstruction,limiter,solver
    logical, intent(in) :: chemistry,viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(inout) :: boundary_change(:),maximum_correction,minimum_theta
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: chemistry_integrator
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: difmag
    character(len=:), allocatable :: transport_limiter
    type(amr_tree_state_3d) :: candidate,temps
    real(dp) :: budget(size(boundary_change)),stage_budget(size(boundary_change)),largest,stage_max,theta,stage_theta
    logical, optional, intent(in) :: subcycling
    logical :: valid,transport_enabled,recursive_steps
    integer :: rank,status,flags(6),low(6),high(6)
    ok=.false.
    call agree_mpi_constant_transport_coefficients(constant_coefficients,comm,constant_agreed,species=species)
    if(.not.constant_agreed) return
    valid=ieee_is_finite(dt).and.dt>0.and..not.(barodiffusion.and..not.diffusion)
    call agree(valid)
    if (.not.valid) return
    recursive_steps=.false.
    if (present(subcycling)) recursive_steps=subcycling
    flags=[merge(1,0,chemistry),merge(1,0,viscosity),merge(1,0,conduction),merge(1,0,diffusion), &
      merge(1,0,barodiffusion),merge(1,0,recursive_steps)]
    call MPI_Allreduce(flags,low,6,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,6,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(low/=high)) return
    if (reconstruction=='characteristic_plm'.or.reconstruction=='characteristic_ppm') transport_limiter=limiter
    call clone_mpi_amr_tree_inputs_3d(tree,owners,size(species),state,temperature,candidate,temps,comm,rank,valid)
    if (.not.valid) return
    budget=0;largest=0;theta=1
    transport_enabled=viscosity.or.conduction.or.diffusion
    if (chemistry) then
      call reaction_half(valid)
      if (.not.valid) return
    end if
    if (transport_enabled) then
      call transport_half(valid)
      if (.not.valid) return
    end if
    stage_budget=0;stage_max=0
    if (recursive_steps) then
      call advance_mpi_amr_tree_subcycle_hydro_3d(tree,owners,species,candidate,temps,spacing,dt,faces, &
        reconstruction,limiter,solver,stage_budget,stage_max,comm,valid,inflow_state,inflow_temperature, &
        outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure,difmag=difmag)
    else
      call advance_mpi_amr_tree_hydro_3d(tree,owners,species,candidate,temps,spacing,dt,faces,reconstruction,limiter,solver, &
        stage_budget,stage_max,comm,valid,inflow_state,inflow_temperature, &
        outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure,difmag=difmag)
    end if
    if (.not.valid) return
    budget=budget+stage_budget;largest=max(largest,stage_max)
    if (transport_enabled) then
      call transport_half(valid)
      if (.not.valid) return
    end if
    if (chemistry) then
      call reaction_half(valid)
      if (.not.valid) return
    end if
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(temps%fields,temperature%fields)
    boundary_change=budget;maximum_correction=largest;minimum_theta=theta;ok=.true.
  contains
    subroutine reaction_half(success)
      logical, intent(out) :: success
      call advance_mpi_amr_tree_chemistry_3d(tree,owners,species,reactions,candidate,temps,.5_dp*dt,rtol,atol, &
        comm,success,chemistry_integrator)
    end subroutine
    subroutine transport_half(success)
      logical, intent(out) :: success
      stage_budget=0;stage_max=0;stage_theta=1
      if (recursive_steps) then
        call advance_mpi_amr_tree_subcycle_transport_3d(tree,owners,species,transport,candidate,temps,spacing, &
          .5_dp*dt,faces,viscosity,conduction,diffusion,barodiffusion,stage_budget,stage_max,stage_theta,comm,success, &
          inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter=transport_limiter, &
          outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure, &
          constant_coefficients=constant_coefficients)
      else
        call advance_mpi_amr_tree_transport_3d(tree,owners,species,transport,candidate,temps,spacing,.5_dp*dt,faces, &
          viscosity,conduction,diffusion,barodiffusion,stage_budget,stage_max,stage_theta,comm,success, &
          inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter=transport_limiter, &
          outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure, &
          constant_coefficients=constant_coefficients)
      end if
      if (.not.success) return
      budget=budget+stage_budget;largest=max(largest,stage_max);theta=min(theta,stage_theta)
    end subroutine
    subroutine agree(condition)
      logical, intent(inout) :: condition
      logical :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      condition=.false.
      if (ierr==MPI_SUCCESS) condition=result
    end subroutine
  end subroutine advance_mpi_amr_tree_split_3d
end module mpi_amr_patch_tree_split_3d_mod
