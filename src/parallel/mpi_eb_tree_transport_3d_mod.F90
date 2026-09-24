module mpi_eb_tree_transport_3d_mod
  use mpi_f08
  use mpi_amr_tree_work_mod, only: tree_work_start,tree_work_stop
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use gas_transport_mod, only: gas_transport_species
  use eb_reactive_transport_3d_mod, only: advance_reactive_eb_transport_patch_euler_3d
  use mpi_eb_tree_transfer_3d_mod, only: average_down_mpi_eb_tree_3d
  use mpi_eb_tree_halo_3d_mod, only: fill_mpi_eb_tree_halo_3d
  use mpi_eb_tree_flux_3d_mod, only: synchronize_mpi_eb_tree_flux_3d
  implicit none
  private
  public :: advance_mpi_eb_tree_transport_3d
contains
  ! Common molecular-transport interval with closed physical diffusion faces.
  subroutine advance_mpi_eb_tree_transport_3d(species,transport,tree,owners,geometry,dt,viscosity,conduction, &
      diffusion,barodiffusion,state,temperature,boundary_change,minimum_theta,comm,ok, &
      inflow_mach,inflow_temperature,inflow_pressure)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(inout) :: boundary_change(:),minimum_theta
    type(MPI_Comm), intent(in) :: comm
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: snapshot,st,candidate,ct,ghost,gt
    type(eb_integrated_flux_3d), allocatable :: flux(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp) :: work_stamp
    real(dp) :: boundary(size(boundary_change)),theta,local_theta,stage_theta,low,high
    integer :: nv,rank,status,node,n(3),c,ios,depth,flags(4),minflags(4),maxflags(4)
    logical :: faces(6)
    logical :: valid,ready
    ok=.false.;nv=reactive_nvar(size(species))
    valid=ieee_is_finite(dt)
    if (valid) valid=dt>0._dp
    valid=valid.and.size(boundary_change)==nv
    call agree(valid,ready)
    if (.not.ready) return
    call MPI_Allreduce(dt,low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(dt,high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS.or.low/=high) return
    flags=merge(1,0,[viscosity,conduction,diffusion,barodiffusion])
    call MPI_Allreduce(flags,minflags,4,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,maxflags,4,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS.or.any(minflags/=maxflags)) return
    snapshot=state;st=temperature
    call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,snapshot,st,comm,ready)
    if (.not.ready) return
    call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,snapshot,st,ghost,gt,comm,ready, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid,ready)
    if (.not.ready) return
    candidate=snapshot;ct=st
    allocate(flux(size(owners)),stat=ios)
    valid=ios==0
    call agree(valid,ready)
    if (.not.ready) return
    local_theta=1._dp
    do node=1,size(owners)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape;depth=tree%nodes(node)%level+1
      faces=[origin(1,node)>0,origin(1,node)+n(1)<domain(1,depth), &
        origin(2,node)>0,origin(2,node)+n(2)<domain(2,depth), &
        origin(3,node)>0,origin(3,node)+n(3)<domain(3,depth)]
      work_stamp=tree_work_start()
      if (present(inflow_mach)) faces=faces.or.(inflow_mach>=0)
      call advance_reactive_eb_transport_patch_euler_3d(species,transport,geometry(node),dt, &
        viscosity,conduction,diffusion,barodiffusion,candidate%fields(node)%values,ct%fields(node)%values(1,:,:,:), &
        flux(node),stage_theta,valid,ghost%fields(node)%values,gt%fields(node)%values(1,:,:,:),faces)
      call tree_work_stop(node,work_stamp)
      if (.not.valid) exit
      local_theta=min(local_theta,stage_theta)
    end do
    call agree(valid.and.ieee_is_finite(local_theta),ready)
    if (.not.ready) return
    call synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,candidate,ct,flux,comm,ready)
    if (.not.ready) return
    boundary=0._dp
    if (rank==owners(1)) then
      n=tree%nodes(1)%shape
      associate(g=>geometry(1),f=>flux(1))
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
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(local_theta,theta,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    call agree(status==MPI_SUCCESS.and.all(ieee_is_finite(boundary)).and.ieee_is_finite(theta),ready)
    if (.not.ready) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(ct%fields,temperature%fields)
    boundary_change=boundary;minimum_theta=theta;ok=.true.
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) result=.false.
    end subroutine
  end subroutine
end module mpi_eb_tree_transport_3d_mod
