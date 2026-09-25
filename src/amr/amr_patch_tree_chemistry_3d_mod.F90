module amr_patch_tree_chemistry_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use reactive_1d_mod, only: reactive_nvar
  use reactive_3d_mod, only: advance_reactive_chemistry_3d,recover_reactive_temperatures_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,average_down_amr_patch_tree_3d
  use amr_patch_tree_hydro_3d_mod, only: advance_amr_patch_tree_hydro_3d
  implicit none
  private
  public :: advance_amr_patch_tree_chemistry_3d,advance_amr_patch_tree_strang_3d
contains
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
