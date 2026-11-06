module amr_patch_tree_chemistry_3d_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use constant_volume_reactor_mod, only: forced_reactor_integrator
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use reactive_1d_mod, only: reactive_nvar
  use rotating_frame_3d_mod, only: rotating_frame_3d,rotating_potential_3d
  use reactive_3d_mod, only: advance_reactive_chemistry_3d,recover_reactive_temperatures_3d,react_sdc_state_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,average_down_amr_patch_tree_3d, &
    amr_tree_integer_layout_3d
  use amr_patch_tree_hydro_3d_mod, only: advance_amr_patch_tree_hydro_3d
  implicit none
  private
  public :: advance_amr_patch_tree_chemistry_3d,advance_amr_patch_tree_strang_3d
  public :: react_sdc_amr_patch_tree_3d
  public :: initialize_amr_tree_reaction_source_3d
contains
  ! Seed the first MOL predictor without advancing the public gas/temperature.
  ! Pinned post_init processes the root, restricts the then-zero fine history,
  ! and subsequently initializes each finer level. Keep that ordering explicit.
  subroutine initialize_amr_tree_reaction_source_3d(tree,species,reactions,state,temperature,dt,subcycling, &
      rtol,atol,reaction_source,ok,evaluate_reaction,internal_energy,frame,domain_lower,root_spacing)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: dt,rtol,atol
    logical, intent(in) :: subcycling
    type(amr_tree_state_3d), intent(inout) :: reaction_source
    logical, intent(out) :: ok
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3),root_spacing(3)
    type(amr_tree_state_3d) :: candidate,temps,rates,sources
    type(amr_tree_state_3d), allocatable :: held
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer :: node,depth,nv
    real(dp) :: level_dt
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    if(.not.state%is_valid(tree).or..not.temperature%is_valid(tree)) return
    if(.not.reaction_source%is_valid(tree)) return
    if(size(state%fields(1)%values,1)/=nv.or.size(temperature%fields(1)%values,1)/=1) return
    if(size(reaction_source%fields(1)%values,1)/=nv-4) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    if(.not.valid) return
    candidate=state;temps=temperature;rates=state;sources=reaction_source
    if(present(internal_energy)) then
      allocate(held);held=internal_energy
    endif
    do node=1,size(tree%nodes)
      rates%fields(node)%values=0;sources%fields(node)%values=0
    enddo
    level_dt=dt
    do depth=0,tree%level_count()-1
      call react_sdc_amr_patch_tree_3d(tree,species,reactions,state,temperature,rates,level_dt,rtol,atol, &
        candidate,temps,sources,valid,target_depth=depth,evaluate_reaction=evaluate_reaction, &
        internal_energy=held,old_internal_energy=internal_energy,frame=frame, &
        domain_lower=domain_lower,root_spacing=root_spacing)
      if(.not.valid) return
      if(depth==0) then
        call average_down_amr_patch_tree_3d(tree,sources,valid)
        if(.not.valid) return
      endif
      if(subcycling.and.depth<tree%level_count()-1) level_dt=level_dt/real(ratios(depth+2),dp)
    enddo
    reaction_source=sources;ok=.true.
  end subroutine initialize_amr_tree_reaction_source_3d

  subroutine react_sdc_amr_patch_tree_3d(tree,species,reactions,old_state,old_temperature,nonreacting_rates, &
      dt,rtol,atol,state,temperature,reaction_source,ok,target_depth,evaluate_reaction,heat_release,internal_energy, &
      frame,domain_lower,root_spacing,old_internal_energy)
    type(amr_tree_state_3d), optional, intent(in) :: old_internal_energy
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3),root_spacing(3)
    type(amr_tree_state_3d), optional, intent(inout) :: heat_release
    type(amr_tree_state_3d), optional, intent(inout) :: internal_energy
    type(amr_tree_state_3d), allocatable :: heat,auxiliary
    procedure(forced_reactor_integrator), optional :: evaluate_reaction
    type(amr_patch_tree_3d), intent(in) :: tree
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_tree_state_3d), intent(in) :: old_state,old_temperature,nonreacting_rates
    real(dp), intent(in) :: dt,rtol,atol
    type(amr_tree_state_3d), intent(inout) :: state,temperature,reaction_source
    logical, intent(out) :: ok
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d) :: candidate,temps,sources
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: frame_h(:,:),node_potential(:,:,:)
    real(dp) :: position(3)
    integer :: node,n(3),nv,depth,i,j,k
    logical :: valid
    real(dp), allocatable :: node_heat(:,:,:),node_auxiliary(:,:,:),node_old_auxiliary(:,:,:)
    ok=.false.;nv=reactive_nvar(size(species))
    if(.not.old_state%is_valid(tree).or..not.old_temperature%is_valid(tree)) return
    if(.not.nonreacting_rates%is_valid(tree)) return
    if(size(old_state%fields(1)%values,1)/=nv.or.size(nonreacting_rates%fields(1)%values,1)/=nv) return
    if(size(old_temperature%fields(1)%values,1)/=1) return
    if(present(old_internal_energy)) then
      if(.not.present(internal_energy)) return
      if(.not.old_internal_energy%is_valid(tree)) return
      if(size(old_internal_energy%fields(1)%values,1)/=1) return
    endif
    if(.not.all(ieee_is_finite([dt,rtol,atol])).or.min(dt,rtol,atol)<=0) return
    if(present(frame).neqv.present(domain_lower)) return
    if(present(frame)) then
      if(.not.present(root_spacing)) return
      if(.not.frame%is_valid()) return
      if(.not.all(ieee_is_finite(domain_lower)).or..not.all(ieee_is_finite(root_spacing))) return
      if(any(root_spacing<=0)) return
      call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
      if(.not.valid) return
      allocate(frame_h(3,size(tree%nodes)))
      frame_h(:,1)=root_spacing
      do node=2,size(tree%nodes)
        frame_h(:,node)=frame_h(:,tree%nodes(node)%parent)/real(tree%nodes(node)%patch%refinement_ratio,dp)
      enddo
    endif
    depth=-1
    if(present(target_depth)) then
      depth=target_depth
      if(depth<0.or.depth>=tree%level_count()) return
      ! Other levels may already be at different times: retain their outputs.
      if(.not.state%is_valid(tree).or..not.temperature%is_valid(tree).or..not.reaction_source%is_valid(tree)) return
      if(size(state%fields(1)%values,1)/=nv.or.size(temperature%fields(1)%values,1)/=1) return
      if(size(reaction_source%fields(1)%values,1)/=nv-4) return
      candidate=state;temps=temperature;sources=reaction_source
    else
      candidate=old_state;temps=old_temperature
      allocate(sources%fields(size(tree%nodes)))
    end if
    if(present(heat_release)) then
      if(.not.heat_release%is_valid(tree)) return
      if(size(heat_release%fields(1)%values,1)/=1) return
      allocate(heat);heat=heat_release
    endif
    if(present(internal_energy)) then
      if(.not.internal_energy%is_valid(tree)) return
      if(size(internal_energy%fields(1)%values,1)/=1) return
      allocate(auxiliary);auxiliary=internal_energy
    endif
    ! SDC evolves each level's valid cells, including covered parents used
    ! during iteration. The hierarchy driver restricts only after reflux;
    ! premature averaging here would change the reaction source estimates.
    do node=1,size(tree%nodes)
      if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
      n=tree%nodes(node)%shape
      if(.not.allocated(sources%fields(node)%values)) allocate(sources%fields(node)%values(nv-4,n(1),n(2),n(3)))
      sources%fields(node)%values=0
      if(present(heat_release)) allocate(node_heat(n(1),n(2),n(3)),source=0._dp)
      if(present(internal_energy)) allocate(node_auxiliary(n(1),n(2),n(3)),source=0._dp)
      if(present(old_internal_energy)) &
        allocate(node_old_auxiliary,source=old_internal_energy%fields(node)%values(1,:,:,:))
      if(present(frame)) then
        allocate(node_potential(n(1),n(2),n(3)))
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          position=domain_lower+(real(origin(:,node),dp)+real([i,j,k],dp)-.5_dp)*frame_h(:,node)
          call rotating_potential_3d(frame,position,node_potential(i,j,k),valid)
          if(.not.valid) return
        enddo
        enddo
        enddo
      endif
      call react_sdc_state_3d(species,reactions,old_state%fields(node)%values, &
        old_temperature%fields(node)%values(1,:,:,:),nonreacting_rates%fields(node)%values,dt,rtol,atol, &
        candidate%fields(node)%values,temps%fields(node)%values(1,:,:,:),sources%fields(node)%values,valid,evaluate_reaction, &
          heat_release=node_heat,internal_energy_density=node_auxiliary,rotating_potential=node_potential, &
          old_internal_energy_density=node_old_auxiliary)
      if(.not.valid) return
      if(allocated(node_potential)) deallocate(node_potential)
      if(allocated(node_old_auxiliary)) deallocate(node_old_auxiliary)
      if(present(internal_energy)) then
        auxiliary%fields(node)%values(1,:,:,:)=node_auxiliary
        deallocate(node_auxiliary)
      endif
      if(present(heat_release)) then
        heat%fields(node)%values(1,:,:,:)=node_heat
        deallocate(node_heat)
      endif
    end do
    if(present(internal_energy)) internal_energy=auxiliary
    if(present(heat_release)) heat_release=heat
    state=candidate;temperature=temps;reaction_source=sources;ok=.true.
  end subroutine react_sdc_amr_patch_tree_3d

  subroutine advance_amr_patch_tree_chemistry_3d(tree,species,reactions,state,temperature,interval,rtol,atol, &
      ok,chemistry_integrator)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: interval,rtol,atol
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: chemistry_integrator
    type(amr_tree_state_3d) :: candidate,temps
    logical, allocatable :: active(:,:,:)
    real(dp), allocatable :: recovered(:,:,:)
    integer :: node,child,n(3),lo(3),hi(3)
    logical :: valid
    ok=.false.
    if (.not.state%is_valid(tree).or..not.temperature%is_valid(tree)) return
    if (size(state%fields(1)%values,1)/=reactive_nvar(size(species))) return
    if (size(temperature%fields(1)%values,1)/=1) return
    if (.not.all(ieee_is_finite([interval,rtol,atol]))) return
    if (interval<0.or.min(rtol,atol)<=0) return
    candidate=state;temps=temperature
    do node=1,size(tree%nodes)
      if (any(temps%fields(node)%values<=0)) return
      n=tree%nodes(node)%shape
      allocate(active(n(1),n(2),n(3)));active=.true.
      ! Only exposed cells own chemistry. Covered parents are supplied by
      ! restriction after every leaf and uncovered ancestor has succeeded.
      do child=node+1,size(tree%nodes)
        if (tree%nodes(child)%parent/=node) cycle
        lo=[tree%nodes(child)%patch%coarse_i_lower,tree%nodes(child)%patch%coarse_j_lower, &
          tree%nodes(child)%patch%coarse_k_lower]
        hi=[tree%nodes(child)%patch%coarse_i_upper,tree%nodes(child)%patch%coarse_j_upper, &
          tree%nodes(child)%patch%coarse_k_upper]
        active(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
      end do
      call advance_reactive_chemistry_3d(species,reactions,candidate%fields(node)%values, &
        temps%fields(node)%values(1,:,:,:),n(1),n(2),n(3),interval,rtol,atol,valid,active,chemistry_integrator)
      if (.not.valid) return
      deallocate(active)
    end do
    call average_down_amr_patch_tree_3d(tree,candidate,valid)
    if (.not.valid) return
    do node=1,size(tree%nodes)
      n=tree%nodes(node)%shape
      allocate(recovered(n(1),n(2),n(3)))
      call recover_reactive_temperatures_3d(species,candidate%fields(node)%values,temps%fields(node)%values(1,:,:,:), &
        n(1),n(2),n(3),recovered,valid)
      if (.not.valid) return
      temps%fields(node)%values(1,:,:,:)=recovered
      deallocate(recovered)
    end do
    state=candidate;temperature=temps;ok=.true.
  end subroutine advance_amr_patch_tree_chemistry_3d

  subroutine advance_amr_patch_tree_strang_3d(tree,species,reactions,state,temperature,root_spacing,dt,faces, &
      reconstruction,limiter,solver,rtol,atol,boundary_change,maximum_correction,ok, &
      chemistry_integrator,inflow_state,inflow_temperature,outlet_pressure,difmag)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: root_spacing(3),dt,rtol,atol
    character(len=*), intent(in) :: faces(6),reconstruction,limiter,solver
    real(dp), intent(inout) :: boundary_change(:),maximum_correction
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: chemistry_integrator
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),outlet_pressure(6)
    real(dp), optional, intent(in) :: difmag
    type(amr_tree_state_3d) :: candidate,temps
    real(dp) :: budget(size(boundary_change)),largest
    logical :: valid
    ok=.false.
    if (.not.ieee_is_finite(dt)) return
    if (dt<=0) return
    candidate=state;temps=temperature
    call advance_amr_patch_tree_chemistry_3d(tree,species,reactions,candidate,temps,.5_dp*dt,rtol,atol, &
      valid,chemistry_integrator)
    if (.not.valid) return
    budget=0;largest=0
    call advance_amr_patch_tree_hydro_3d(tree,species,candidate,temps,root_spacing,dt,faces, &
      reconstruction,limiter,solver,budget,largest,valid,inflow_state,inflow_temperature,outlet_pressure,difmag=difmag)
    if (.not.valid) return
    call advance_amr_patch_tree_chemistry_3d(tree,species,reactions,candidate,temps,.5_dp*dt,rtol,atol, &
      valid,chemistry_integrator)
    if (.not.valid) return
    state=candidate;temperature=temps;boundary_change=budget;maximum_correction=largest;ok=.true.
  end subroutine advance_amr_patch_tree_strang_3d
end module amr_patch_tree_chemistry_3d_mod
