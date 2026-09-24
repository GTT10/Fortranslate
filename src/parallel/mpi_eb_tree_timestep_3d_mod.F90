module mpi_eb_tree_timestep_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use reactive_eb_cfl_3d_mod, only: reactive_eb_signal_speeds_3d,compute_reactive_eb_euler_cfl_timestep_3d
  use eb_reactive_transport_3d_mod, only: reactive_eb_transport_patch_timestep_3d,reactive_eb_halo_bounds_3d
  use mpi_eb_tree_halo_3d_mod, only: fill_mpi_eb_tree_halo_3d
  use mpi_eb_tree_transfer_3d_mod, only: average_down_mpi_eb_tree_3d
  implicit none
  private
  public :: mpi_eb_tree_timestep_3d
contains
  ! Root-interval bound including hidden parents and all remote donor maxima.
  ! Only owner-local fields are scanned; geometry remains replicated.
  subroutine mpi_eb_tree_timestep_3d(species,transport,tree,owners,geometry,state,temperature, &
      hydro_cfl,transport_cfl,transport_enabled,viscosity,conduction,diffusion,dt,comm,ok,subcycling, &
      inflow_mach,inflow_temperature,inflow_pressure)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: hydro_cfl,transport_cfl
    logical, intent(in) :: transport_enabled,viscosity,conduction,diffusion
    real(dp), intent(inout) :: dt
    type(MPI_Comm), intent(in) :: comm
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    logical, intent(out) :: ok
    logical, optional, intent(in) :: subcycling
    logical :: recursive_steps
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    type(amr_tree_state_3d) :: snapshot,st,ghost,gt
    real(dp) :: parameters(2),low(2),high(2),speeds(4),maxima(5),donors(5),interval,diffusivity,candidate,result,factor
    integer :: flags(5),minflags(5),maxflags(5),rank,status,node,level,supplied,any_supplied
    logical :: valid,ready,active,any_active
    ok=.false.;parameters=[hydro_cfl,transport_cfl]
    valid=all(ieee_is_finite(parameters))
    if (valid) valid=hydro_cfl>0._dp.and.hydro_cfl<=1._dp.and.transport_cfl>0._dp.and.transport_cfl<=.5_dp
    call agree(valid,ready)
    if (.not.ready) return
    call MPI_Allreduce(parameters,low,2,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(parameters,high,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS.or.any(low/=high)) return
    recursive_steps=.false.
    if (present(subcycling)) recursive_steps=subcycling
    flags=merge(1,0,[transport_enabled,viscosity,conduction,diffusion,recursive_steps])
    call MPI_Allreduce(flags,minflags,5,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,maxflags,5,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS.or.any(minflags/=maxflags)) return
    snapshot=state;st=temperature
    call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,snapshot,st,comm,ready)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    if (recursive_steps) then
      call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
      call agree(valid,ready)
      if (.not.ready) return
    end if
    supplied=merge(1,0,present(inflow_mach).or.present(inflow_temperature).or.present(inflow_pressure))
    call MPI_Allreduce(supplied,any_supplied,1,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any_supplied>0) then
      call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,snapshot,st,ghost,gt,comm,ready, &
        inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure)
      if (.not.ready) return
    end if
    maxima=0._dp;active=.false.
    do node=1,size(owners)
      if (owners(node)/=rank) cycle
      if (.not.any(geometry(node)%volume_fraction>0._dp)) cycle
      active=.true.
      call reactive_eb_signal_speeds_3d(species,snapshot%fields(node)%values,st%fields(node)%values(1,:,:,:), &
        geometry(node),speeds,valid)
      if (.not.valid) exit
      maxima(:4)=max(maxima(:4),speeds)
      if (any_supplied>0) then
        call reactive_eb_halo_bounds_3d(species,transport,geometry(node),ghost%fields(node)%values, &
          gt%fields(node)%values(1,:,:,:),transport_enabled,viscosity,conduction,diffusion,speeds,diffusivity,valid)
        if (.not.valid) exit
        maxima(:4)=max(maxima(:4),speeds);maxima(5)=max(maxima(5),diffusivity)
      end if
      if (transport_enabled) then
        call reactive_eb_transport_patch_timestep_3d(species,transport,snapshot%fields(node)%values, &
          st%fields(node)%values(1,:,:,:),geometry(node),transport_cfl,viscosity,conduction,diffusion, &
          interval,diffusivity,valid)
        if (.not.valid) exit
        maxima(5)=max(maxima(5),diffusivity)
      end if
    end do
    call agree(valid,ready)
    if (.not.ready) return
    call MPI_Allreduce(active,any_active,1,MPI_LOGICAL,MPI_LOR,comm,status)
    if (status/=MPI_SUCCESS.or..not.any_active) return
    call MPI_Allreduce(maxima,donors,5,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    candidate=huge(1._dp)
    do node=1,size(owners)
      if (owners(node)/=rank) cycle
      if (.not.any(geometry(node)%volume_fraction>0._dp)) cycle
      call compute_reactive_eb_euler_cfl_timestep_3d(species,snapshot%fields(node)%values, &
        st%fields(node)%values(1,:,:,:),geometry(node),hydro_cfl,interval,valid,donors(:4))
      if (.not.valid) exit
      factor=1._dp
      if (recursive_steps) then
        level=tree%nodes(node)%level+1
        factor=real(domain(1,level),dp)/real(domain(1,1),dp)
      end if
      if (interval<candidate/factor) candidate=interval*factor
      if (transport_enabled) then
        call reactive_eb_transport_patch_timestep_3d(species,transport,snapshot%fields(node)%values, &
          st%fields(node)%values(1,:,:,:),geometry(node),transport_cfl,viscosity,conduction,diffusion, &
          interval,diffusivity,valid,donors(5))
        if (.not.valid) exit
        ! Retain the serial conservative full-interval bound for split halves.
        factor=factor*factor
        if (interval<candidate/factor) candidate=interval*factor
      end if
    end do
    call agree(valid,ready)
    if (.not.ready) return
    call MPI_Allreduce(candidate,result,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    dt=result;ok=.true.
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) result=.false.
    end subroutine
  end subroutine
end module mpi_eb_tree_timestep_3d_mod
