module les_amr_sdc_mod
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use amr_patch_tree_split_3d_mod, only: advance_amr_patch_tree_sdc_3d
  use les_amr_patch_tree_mod, only: compute_amr_patch_tree_les_rhs_3d
  implicit none
  private
  public :: advance_amr_patch_tree_les_sdc_3d
contains
  subroutine advance_amr_patch_tree_les_sdc_3d(tree,species,reactions,transport,state,temperature,spacing,dt,faces, &
      iterations,chemistry,transport_flags,force,rtol,atol,cs,prt,sct,boundary_change,source_change, &
      maximum_correction,minimum_theta,maximum_diffusivity,ok,reaction_source,model,ci, &
      test_model,test_ratio,grid_ratio,harmonic_mean,inflow_state,inflow_temperature,wall_temperature,wall_velocity, &
      outlet_pressure,inflow_mach,inflow_pressure,difmag)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: spacing(3),dt,force(3),rtol,atol,cs,prt,sct
    character(len=*), intent(in) :: faces(6)
    integer, intent(in) :: iterations
    logical, intent(in) :: chemistry,transport_flags(4)
    real(dp), intent(inout) :: boundary_change(:),source_change(:),maximum_correction,minimum_theta,maximum_diffusivity
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source
    character(len=*), optional, intent(in) :: model,test_model
    real(dp), optional, intent(in) :: ci,difmag,inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6)
    integer, optional, intent(in) :: test_ratio,grid_ratio
    logical, optional, intent(in) :: harmonic_mean
    real(dp) :: largest
    largest=0
    call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,state,temperature,spacing,dt,faces, &
      iterations,chemistry,transport_flags,force,rtol,atol,boundary_change,source_change,maximum_correction,minimum_theta, &
      ok,reaction_source,inflow_state,inflow_temperature,wall_temperature,wall_velocity,outlet_pressure, &
      inflow_mach,inflow_pressure,difmag,evaluate_diffusion=evaluate)
    if(ok) maximum_diffusivity=largest
  contains
    subroutine evaluate(u,t,depth,fractions,rates,flux,good)
      type(amr_tree_state_3d), intent(in) :: u,t
      integer, intent(in) :: depth
      real(dp), optional, intent(in) :: fractions(:)
      type(amr_tree_state_3d), intent(inout) :: rates
      type(amr_tree_flux_3d), allocatable, intent(inout) :: flux(:)
      logical, intent(out) :: good
      real(dp) :: diffusion
      good=.false.
      if(depth/=-1.or.present(fractions)) return
      call compute_amr_patch_tree_les_rhs_3d(tree,species,u,t,spacing,faces,cs,prt,sct,rates,flux,diffusion,good, &
        model,ci,test_model,test_ratio,grid_ratio,harmonic_mean,inflow_state,inflow_temperature, &
        wall_temperature,wall_velocity,outlet_pressure,inflow_mach,inflow_pressure)
      if(good) largest=max(largest,diffusion)
    end subroutine evaluate
  end subroutine advance_amr_patch_tree_les_sdc_3d
end module les_amr_sdc_mod
