module mpi_eb_tree_hydro_3d_mod
  use mpi_f08
  use mpi_amr_tree_work_mod, only: tree_work_start,tree_work_stop
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_reactive_hydro_3d_mod, only: advance_reactive_eb_euler_3d
  use eb_reactive_wall_flux_3d_mod, only: reactive_eb_slip_wall_source_3d
  use mpi_eb_tree_transfer_3d_mod, only: average_down_mpi_eb_tree_3d
  use mpi_eb_tree_halo_3d_mod, only: fill_mpi_eb_tree_halo_3d
  use mpi_eb_tree_flux_3d_mod, only: synchronize_mpi_eb_tree_flux_3d
  implicit none
  private
  public :: advance_mpi_eb_tree_hydro_3d
contains
  ! Common Euler interval, PCM face donors and physical outflow, matching serial.
  subroutine advance_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,solver,dt,state,temperature, &
      boundary_change,wall_change,comm,ok, &
      inflow_mach,inflow_temperature,inflow_pressure)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:)
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(inout) :: boundary_change(:),wall_change(:)
    type(MPI_Comm), intent(in) :: comm
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: snapshot,st,candidate,ct,ghost,gt
    type(eb_integrated_flux_3d), allocatable :: flux(:)
    real(dp), allocatable :: source(:,:,:,:)
    real(dp) :: work_stamp
    real(dp) :: boundary(size(boundary_change)),wall(size(wall_change)),local_wall(size(wall_change)),low,high
    logical, allocatable :: exposed(:,:,:)
    integer :: nv,rank,status,node,child,n(3),lo(3),hi(3),c,ios
    character(len=32) :: reference
    logical :: valid,ready
    ok=.false.;nv=reactive_nvar(size(species))
    valid=ieee_is_finite(dt)
    if (valid) valid=dt>0._dp
    valid=valid.and.size(boundary_change)==nv.and.size(wall_change)==nv.and.len_trim(solver)<=len(reference)
    call agree(valid,ready)
    if (.not.ready) return
    call MPI_Allreduce(dt,low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(dt,high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS.or.low/=high) return
    reference=solver
    call MPI_Bcast(reference,len(reference),MPI_CHARACTER,0,comm,status)
    call agree(status==MPI_SUCCESS.and.trim(solver)==trim(reference),ready)
    if (.not.ready) return
    snapshot=state;st=temperature
    call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,snapshot,st,comm,ready)
    if (.not.ready) return
    call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,snapshot,st,ghost,gt,comm,ready, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure)
    if (.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    candidate=snapshot;ct=st
    allocate(flux(size(owners)),stat=ios)
    valid=ios==0
    call agree(valid,ready)
    if (.not.ready) return
    local_wall=0._dp
    do node=1,size(owners)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      work_stamp=tree_work_start()
      call advance_reactive_eb_euler_3d(species,snapshot%fields(node)%values,st%fields(node)%values(1,:,:,:), &
        geometry(node),solver,dt,candidate%fields(node)%values,ct%fields(node)%values(1,:,:,:), &
        valid,ghost%fields(node)%values,gt%fields(node)%values(1,:,:,:),flux(node))
      call tree_work_stop(node,work_stamp)
      if (.not.valid) exit
      allocate(source(nv,n(1),n(2),n(3)),exposed(n(1),n(2),n(3)),stat=ios)
      valid=ios==0
      if (.not.valid) exit
      exposed=.true.
      do child=node+1,size(owners)
        if (tree%nodes(child)%parent/=node) cycle
        associate(p=>tree%nodes(child)%patch)
        lo=[p%coarse_i_lower,p%coarse_j_lower,p%coarse_k_lower]
        hi=[p%coarse_i_upper,p%coarse_j_upper,p%coarse_k_upper]
        end associate
        exposed(lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))=.false.
      end do
      call reactive_eb_slip_wall_source_3d(species,snapshot%fields(node)%values,st%fields(node)%values(1,:,:,:), &
        geometry(node),source,valid)
      if (.not.valid) exit
      associate(g=>geometry(node))
      do c=1,nv
        local_wall(c)=local_wall(c)+dt*g%dx*g%dy*g%dz*sum(source(c,:,:,:)*g%volume_fraction,mask=exposed)
      end do
      end associate
      deallocate(source,exposed)
    end do
    call agree(valid.and.all(ieee_is_finite(local_wall)),ready)
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
    call MPI_Allreduce(local_wall,wall,nv,MPI_DOUBLE_PRECISION,MPI_SUM,comm,status)
    call agree(status==MPI_SUCCESS.and.all(ieee_is_finite(boundary)).and.all(ieee_is_finite(wall)),ready)
    if (.not.ready) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(ct%fields,temperature%fields)
    boundary_change=boundary;wall_change=wall;ok=.true.
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if (ierr/=MPI_SUCCESS) result=.false.
    end subroutine
  end subroutine
end module mpi_eb_tree_hydro_3d_mod
