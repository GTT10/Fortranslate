module les_mpi_eb_sdc_mod
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,eb_regular_cell_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use mpi_eb_tree_transfer_3d_mod, only: validate_mpi_eb_tree_geometry_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use les_mpi_amr_tree_mod, only: compute_mpi_amr_tree_les_rhs_3d
  implicit none
  private
  public :: compute_mpi_regular_eb_les_sdc_3d
contains
  ! The pinned LES model applies only to a completely regular geometry.
  ! Convert separate EB temporal donors and zero-based integrated face fluxes.
  subroutine compute_mpi_regular_eb_les_sdc_3d(tree,owners,geometry,species,state,temperature,spacing,faces, &
      interval,depth,cs,prt,sct,rates,flux,maximum_diffusivity,comm,ok, &
      ancestor_start,ancestor_end,temperature_start,temperature_end,fractions,model,ci, &
      test_model,test_ratio,grid_ratio,harmonic_mean,inflow_state,inflow_mach,inflow_temperature,inflow_pressure, &
      wall_temperature,wall_velocity)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:),depth
    type(eb_geometry_3d), intent(in) :: geometry(:)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: spacing(3),interval,cs,prt,sct
    character(len=*), intent(in) :: faces(6)
    type(amr_tree_state_3d), intent(inout) :: rates
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    real(dp), intent(inout) :: maximum_diffusivity
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: fractions(:)
    character(len=*), optional, intent(in) :: model,test_model
    real(dp), optional, intent(in) :: ci,inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: wall_temperature(6),wall_velocity(:,:)
    integer, optional, intent(in) :: test_ratio,grid_ratio
    logical, optional, intent(in) :: harmonic_mean
    type(amr_tree_state_3d), allocatable :: packed_start,packed_end
    type(amr_tree_state_3d) :: candidate
    type(amr_tree_flux_3d), allocatable :: density_flux(:)
    type(eb_integrated_flux_3d), allocatable :: result(:)
    real(dp), allocatable :: stage_fractions(:)
    real(dp) :: largest,low_dt,high_dt
    integer :: rank,status,node,n(3),nv,ios,controls(6),low(6),high(6)
    logical :: valid,ready
    ok=.false.;nv=5+size(species)
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    controls=[depth,merge(1,0,present(ancestor_start)),merge(1,0,present(ancestor_end)), &
      merge(1,0,present(temperature_start)),merge(1,0,present(temperature_end)),merge(1,0,present(fractions))]
    call MPI_Allreduce(controls,low,6,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,6,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(low/=high).or.depth< -1.or.depth>=tree%level_count()) return
    valid=ieee_is_finite(interval).and.interval>0
    do node=1,size(geometry)
      valid=valid.and.all(geometry(node)%cell_type==eb_regular_cell_3d)
      valid=valid.and.all(geometry(node)%volume_fraction==1)
      valid=valid.and.all(geometry(node)%x_face_fraction==1).and.all(geometry(node)%y_face_fraction==1).and. &
        all(geometry(node)%z_face_fraction==1)
    end do
    call agree(valid)
    if(.not.valid) return
    call MPI_Allreduce(interval,low_dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(interval,high_dt,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(low_dt/=high_dt) return
    if(depth>0) then
      if(any(controls(2:)/=1)) return
      valid=size(fractions)>=depth
      call agree(valid)
      if(.not.valid) return
      call validate_mpi_amr_tree_local_3d(tree,owners,ancestor_start,nv,comm,ready,maximum_depth=depth-1)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,ancestor_end,nv,comm,ready,maximum_depth=depth-1)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,temperature_start,1,comm,ready,maximum_depth=depth-1)
      if(.not.ready) return
      call validate_mpi_amr_tree_local_3d(tree,owners,temperature_end,1,comm,ready,maximum_depth=depth-1)
      if(.not.ready) return
      allocate(packed_start,packed_end,stat=ios)
      valid=ios==0;call agree(valid)
      if(.not.valid) return
      allocate(packed_start%fields(size(owners)),packed_end%fields(size(owners)),stage_fractions(depth),stat=ios)
      valid=ios==0;call agree(valid)
      if(.not.valid) return
      stage_fractions=fractions(:depth)
      do node=1,size(owners)
        if(owners(node)/=rank.or.tree%nodes(node)%level>=depth) cycle
        n=tree%nodes(node)%shape
        allocate(packed_start%fields(node)%values(nv+1,n(1),n(2),n(3)), &
          packed_end%fields(node)%values(nv+1,n(1),n(2),n(3)),stat=ios)
        if(ios/=0) exit
        packed_start%fields(node)%values(:nv,:,:,:)=ancestor_start%fields(node)%values
        packed_start%fields(node)%values(nv+1,:,:,:)=temperature_start%fields(node)%values(1,:,:,:)
        packed_end%fields(node)%values(:nv,:,:,:)=ancestor_end%fields(node)%values
        packed_end%fields(node)%values(nv+1,:,:,:)=temperature_end%fields(node)%values(1,:,:,:)
      end do
      valid=ios==0;call agree(valid)
      if(.not.valid) return
    end if
    call compute_mpi_amr_tree_les_rhs_3d(tree,owners,species,state,temperature,spacing,faces,cs,prt,sct, &
      candidate,density_flux,largest,comm,ready,inflow_state=inflow_state,inflow_temperature=inflow_temperature, &
      wall_temperature=wall_temperature,wall_velocity=wall_velocity,model=model,ci=ci,test_model=test_model, &
      test_ratio=test_ratio,grid_ratio=grid_ratio,harmonic_mean=harmonic_mean,inflow_mach=inflow_mach, &
      inflow_pressure=inflow_pressure,target_depth=depth,coarse_start=packed_start,coarse_end=packed_end,fractions=stage_fractions)
    if(.not.ready) return
    allocate(result(size(owners)),stat=ios)
    valid=ios==0;call agree(valid)
    if(.not.valid) return
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
      n=tree%nodes(node)%shape
      allocate(result(node)%x(nv,0:n(1),n(2),n(3)),result(node)%y(nv,n(1),0:n(2),n(3)), &
        result(node)%z(nv,n(1),n(2),0:n(3)),stat=ios)
      if(ios/=0) exit
      result(node)%x=interval*density_flux(node)%x
      result(node)%y=interval*density_flux(node)%y
      result(node)%z=interval*density_flux(node)%z
    end do
    valid=ios==0;call agree(valid)
    if(.not.valid) return
    call move_alloc(candidate%fields,rates%fields)
    call move_alloc(result,flux)
    maximum_diffusivity=largest;ok=.true.
  contains
    subroutine agree(condition)
      logical, intent(inout) :: condition
      logical :: global
      integer :: ierr
      call MPI_Allreduce(condition,global,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      condition=.false.
      if(ierr==MPI_SUCCESS) condition=global
    end subroutine
  end subroutine
end module
