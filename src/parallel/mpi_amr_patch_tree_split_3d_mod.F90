module mpi_amr_patch_tree_split_3d_mod
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use mpi_amr_patch_tree_chemistry_3d_mod, only: advance_mpi_amr_tree_chemistry_3d,clone_mpi_amr_tree_inputs_3d
  use mpi_amr_patch_tree_hydro_3d_mod, only: advance_mpi_amr_tree_hydro_3d
  use mpi_amr_patch_tree_transport_3d_mod, only: advance_mpi_amr_tree_transport_3d
  use mpi_amr_patch_tree_subcycle_3d_mod, only: advance_mpi_amr_tree_subcycle_hydro_3d
  use mpi_amr_patch_tree_transport_subcycle_3d_mod, only: advance_mpi_amr_tree_subcycle_transport_3d
  implicit none
  private
  public :: advance_mpi_amr_tree_split_3d
contains
  ! R/2 - T/2 - H - T/2 - R/2, over one coarse interval, optionally with recursive flow/transport steps.
  ! Disabled operators are skipped. All stage outputs remain private until
  ! the complete composition succeeds, including the final chemistry solve.
  subroutine advance_mpi_amr_tree_split_3d(tree,owners,species,reactions,transport,state,temperature,spacing,dt,faces, &
      reconstruction,limiter,solver,chemistry,viscosity,conduction,diffusion,barodiffusion,rtol,atol, &
      boundary_change,maximum_correction,minimum_theta,comm,ok,chemistry_integrator, &
      inflow_state,inflow_temperature,wall_temperature,wall_velocity,subcycling,outlet_pressure,inflow_mach,inflow_pressure)
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
    character(len=:), allocatable :: transport_limiter
    type(amr_tree_state_3d) :: candidate,temps
    real(dp) :: budget(size(boundary_change)),stage_budget(size(boundary_change)),largest,stage_max,theta,stage_theta
    logical, optional, intent(in) :: subcycling
    logical :: valid,transport_enabled,recursive_steps
    integer :: rank,status,flags(6),low(6),high(6)
    ok=.false.
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
        outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
    else
      call advance_mpi_amr_tree_hydro_3d(tree,owners,species,candidate,temps,spacing,dt,faces,reconstruction,limiter,solver, &
        stage_budget,stage_max,comm,valid,inflow_state,inflow_temperature, &
        outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
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
          outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
      else
        call advance_mpi_amr_tree_transport_3d(tree,owners,species,transport,candidate,temps,spacing,.5_dp*dt,faces, &
          viscosity,conduction,diffusion,barodiffusion,stage_budget,stage_max,stage_theta,comm,success, &
          inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter=transport_limiter, &
          outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
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
