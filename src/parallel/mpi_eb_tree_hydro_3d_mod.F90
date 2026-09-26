module mpi_eb_tree_hydro_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use mpi_amr_tree_work_mod, only: tree_work_start,tree_work_stop
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_reactive_hydro_3d_mod, only: advance_reactive_eb_euler_3d, &
    advance_reactive_eb_plm_ctu_3d,compute_reactive_eb_flattening_3d,eb_plm_boundary_3d,compute_sdc_eb_hydro_3d
  use eb_reactive_wall_flux_3d_mod, only: reactive_eb_slip_wall_source_3d
  use mpi_eb_tree_transfer_3d_mod, only: average_down_mpi_eb_tree_3d,validate_mpi_eb_tree_geometry_3d
  use eb_reactive_two_level_3d_mod, only: matching_eb_metrics_3d
  use mpi_eb_tree_halo_3d_mod, only: fill_mpi_eb_tree_halo_3d,fill_sdc_mpi_eb_tree_source_halo_3d
  use mpi_eb_tree_flux_3d_mod, only: synchronize_mpi_eb_tree_flux_3d
  implicit none
  private
  public :: advance_mpi_eb_tree_hydro_3d,compute_sdc_mpi_eb_tree_hydro_3d
contains
  ! Common interval with PCM or an explicit six-layer PLM CTU geometry cache.
  subroutine advance_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,solver,dt,state,temperature, &
      boundary_change,wall_change,comm,ok, &
      inflow_mach,inflow_temperature,inflow_pressure,difmag,extended_geometry,boundary_faces,inflow_state)
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
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: inflow_state(:,:)
    real(dp), optional, intent(in) :: difmag
    type(eb_geometry_3d), optional, intent(in) :: extended_geometry(:)
    type(eb_geometry_3d) :: center,working
    type(eb_plm_boundary_3d) :: physical_boundary
    type(amr_tree_state_3d) :: metrics
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: flat(:,:,:)
    logical :: ctu,reference_ctu
    real(dp) :: coefficient
    logical, intent(out) :: ok
    type(amr_tree_state_3d) :: snapshot,st,candidate,ct,ghost,gt,vf
    type(eb_integrated_flux_3d), allocatable :: flux(:)
    real(dp), allocatable :: source(:,:,:,:)
    real(dp) :: work_stamp
    real(dp) :: boundary(size(boundary_change)),wall(size(wall_change)),local_wall(size(wall_change)),low,high
    logical, allocatable :: exposed(:,:,:)
    integer :: nv,rank,status,node,child,n(3),lo(3),hi(3),c,ios
    character(len=32) :: reference
    logical :: valid,ready
    ok=.false.;nv=reactive_nvar(size(species))
    coefficient=.1_dp
    if(present(difmag)) coefficient=difmag
    valid=ieee_is_finite(dt).and.ieee_is_finite(coefficient).and.coefficient>=0
    if (valid) valid=dt>0._dp
    valid=valid.and.size(boundary_change)==nv.and.size(wall_change)==nv.and.len_trim(solver)<=len(reference)
    call agree(valid,ready)
    if (.not.ready) return
    call MPI_Allreduce(dt,low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(dt,high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS.or.low/=high) return
    call MPI_Allreduce(coefficient,low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(coefficient,high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.low/=high) return
    reference=solver
    call MPI_Bcast(reference,len(reference),MPI_CHARACTER,0,comm,status)
    call agree(status==MPI_SUCCESS.and.trim(solver)==trim(reference),ready)
    if (.not.ready) return
    ctu=present(extended_geometry);reference_ctu=ctu
    call MPI_Bcast(reference_ctu,1,MPI_LOGICAL,0,comm,status)
    call agree(status==MPI_SUCCESS.and.(ctu.eqv.reference_ctu),ready)
    if(.not.ready) return
    if(ctu) then
      call agree(solver=='pelec',ready)
      if(.not.ready) return
      call validate_mpi_eb_tree_geometry_3d(tree,owners,extended_geometry,nv,comm,ready,padding=6)
      if(.not.ready) return
    end if
    snapshot=state;st=temperature
    call average_down_mpi_eb_tree_3d(species,tree,owners,geometry,snapshot,st,comm,ready)
    if (.not.ready) return
    if(ctu) then
      call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
      if(valid) valid=all(domain<=int(huge(1)-4,int64))
      call agree(valid,ready)
      if(.not.ready) return
      call MPI_Comm_rank(comm,rank,status)
      if(status/=MPI_SUCCESS) return
      allocate(metrics%fields(size(owners)))
      do node=1,size(owners)
        n=tree%nodes(node)%shape
        call extract_eb_geometry_3d(extended_geometry(node),[7,7,7],n+6,center,valid)
        if(valid) valid=matching_eb_metrics_3d(center,geometry(node))
        if(.not.valid) exit
        if(owners(node)/=rank) cycle
        allocate(metrics%fields(node)%values(1,-5:n(1)+6,-5:n(2)+6,-5:n(3)+6))
        metrics%fields(node)%values(1,:,:,:)=extended_geometry(node)%volume_fraction
      end do
      call agree(valid,ready)
      if(.not.ready) return
      call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,snapshot,st,ghost,gt,comm,ready, &
        inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state, &
        ghost_volume_fraction=vf,halo_width=6,target_volume_fraction=metrics)
    else if(coefficient>0) then
    call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,snapshot,st,ghost,gt,comm,ready, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state,ghost_volume_fraction=vf)
    else
    call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,snapshot,st,ghost,gt,comm,ready, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        boundary_faces=boundary_faces,inflow_state=inflow_state)
    end if
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
      if(ctu) then
        call extract_eb_geometry_3d(extended_geometry(node),[4,4,4],n+9,working,valid)
        if(.not.valid) exit
        allocate(flat(0:n(1)+7,0:n(2)+7,0:n(3)+7))
        call compute_reactive_eb_flattening_3d(species,ghost%fields(node)%values,gt%fields(node)%values(1,:,:,:), &
          vf%fields(node)%values(1,:,:,:),n+6,flat,valid)
        if(.not.valid) exit
        physical_boundary%lower=4-int(origin(:,node))
        physical_boundary%upper=3+int(domain(:,tree%nodes(node)%level+1)-origin(:,node))
        physical_boundary%faces='outflow'
        if(present(boundary_faces)) physical_boundary%faces=boundary_faces
        if(present(inflow_mach)) then
          where(inflow_mach>=0) physical_boundary%faces='inflow'
        end if
        call advance_reactive_eb_plm_ctu_3d(species,working,ghost%fields(node)%values, &
          gt%fields(node)%values(1,:,:,:),vf%fields(node)%values(1,:,:,:),dt,flat,candidate%fields(node)%values, &
          ct%fields(node)%values(1,:,:,:),valid,difmag=coefficient,integrated_flux=flux(node),boundary=physical_boundary)
        if(.not.valid) exit
        do c=1,nv
          where(geometry(node)%volume_fraction<=0) candidate%fields(node)%values(c,:,:,:)=snapshot%fields(node)%values(c,:,:,:)
        end do
        where(geometry(node)%volume_fraction<=0) ct%fields(node)%values(1,:,:,:)=st%fields(node)%values(1,:,:,:)
        deallocate(flat)
      else if(coefficient>0) then
      call advance_reactive_eb_euler_3d(species,snapshot%fields(node)%values,st%fields(node)%values(1,:,:,:), &
        geometry(node),solver,dt,candidate%fields(node)%values,ct%fields(node)%values(1,:,:,:), &
        valid,ghost%fields(node)%values,gt%fields(node)%values(1,:,:,:),flux(node), &
        difmag=coefficient,ghost_volume_fraction=vf%fields(node)%values(1,:,:,:))
      else
      call advance_reactive_eb_euler_3d(species,snapshot%fields(node)%values,st%fields(node)%values(1,:,:,:), &
        geometry(node),solver,dt,candidate%fields(node)%values,ct%fields(node)%values(1,:,:,:), &
        valid,ghost%fields(node)%values,gt%fields(node)%values(1,:,:,:),flux(node),difmag=coefficient)
      end if
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
    call synchronize_mpi_eb_tree_flux_3d(species,tree,owners,geometry,candidate,ct,flux,comm,ready, &
      boundary_faces=boundary_faces)
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
  ! SDC CTU stage: owner-local rates and integrated fluxes, without reflux or restriction.
  subroutine compute_sdc_mpi_eb_tree_hydro_3d(species,tree,owners,geometry,state,temperature,predictor, &
      extended_geometry,dt,hydro_rate,flux,comm,ok,difmag,boundary_faces,inflow_state, &
      inflow_mach,inflow_temperature,inflow_pressure, &
      target_depth,ancestor_start,ancestor_end,temperature_start,temperature_end,ancestor_fraction, &
      source_start,source_end,source_fraction,predictor_halo)
    type(nasa7_species), intent(in) :: species(:)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(eb_geometry_3d), intent(in) :: geometry(:),extended_geometry(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature,predictor
    real(dp), intent(in) :: dt
    type(amr_tree_state_3d), intent(inout) :: hydro_rate
    type(eb_integrated_flux_3d), allocatable, intent(inout) :: flux(:)
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d), optional, intent(in) :: ancestor_start,ancestor_end,temperature_start,temperature_end
    real(dp), optional, intent(in) :: ancestor_fraction(:)
    integer :: selected_depth,min_depth,max_depth
    type(amr_tree_state_3d), optional, intent(in) :: source_start,source_end,predictor_halo
    real(dp), optional, intent(in) :: source_fraction(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: difmag,inflow_state(:,:),inflow_mach(6),inflow_temperature(6),inflow_pressure(6)
    character(len=*), optional, intent(in) :: boundary_faces(6)
    type(amr_tree_state_3d) :: rates,metrics,source_metrics,ghost,gt,vf,source
    type(eb_integrated_flux_3d), allocatable :: candidate_flux(:)
    type(eb_geometry_3d) :: center,working
    type(eb_plm_boundary_3d) :: physical_boundary
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: flat(:,:,:)
    real(dp) :: coefficient,controls(2),low(2),high(2),work_stamp
    real(dp) :: boundary(reactive_nvar(size(species))),wall(reactive_nvar(size(species)))
    integer :: nv,rank,status,ios,node,n(3),flag,loflag,hiflag,i,j,k,offset(3)
    character(len=32) :: faces(6)
    logical :: valid,ready
    ok=.false.;nv=reactive_nvar(size(species))
    coefficient=.1_dp
    if(present(difmag)) coefficient=difmag
    controls=[dt,coefficient]
    call agree(all(ieee_is_finite(controls)).and.dt>0.and.coefficient>=0,ready)
    if(.not.ready) return
    call MPI_Allreduce(controls,low,2,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    flag=merge(1,0,present(predictor_halo))
    call MPI_Allreduce(flag,loflag,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flag,hiflag,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.loflag/=hiflag) return
    selected_depth=-1
    if(present(target_depth)) selected_depth=target_depth
    valid=selected_depth>=-1
    if(present(target_depth)) valid=selected_depth>=0
    if(present(ancestor_start).or.present(ancestor_end).or.present(temperature_start).or. &
      present(temperature_end).or.present(ancestor_fraction)) valid=valid.and.selected_depth>=0
    if(present(source_start).or.present(source_end).or.present(source_fraction)) valid=valid.and.selected_depth>=0
    call agree(valid,ready)
    if(.not.ready) return
    call MPI_Allreduce(selected_depth,min_depth,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(selected_depth,max_depth,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.min_depth/=max_depth) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,geometry,nv,comm,ready)
    if(.not.ready) return
    call validate_mpi_eb_tree_geometry_3d(tree,owners,extended_geometry,nv,comm,ready,padding=6)
    if(.not.ready) return
    call agree(selected_depth<tree%level_count(),ready)
    if(.not.ready) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    if(valid) valid=all(domain<=int(huge(1)-4,int64))
    call agree(valid,ready)
    if(.not.ready) return
    call MPI_Comm_rank(comm,rank,status)
    if(status/=MPI_SUCCESS) return
    allocate(metrics%fields(size(owners)),source_metrics%fields(size(owners)), &
      rates%fields(size(owners)),candidate_flux(size(owners)),stat=ios)
    call agree(ios==0,ready)
    if(.not.ready) return
    do node=1,size(owners)
      n=tree%nodes(node)%shape
      call extract_eb_geometry_3d(extended_geometry(node),[7,7,7],n+6,center,valid)
      if(valid) valid=matching_eb_metrics_3d(center,geometry(node))
      if(.not.valid) exit
      if(owners(node)/=rank) cycle
      allocate(metrics%fields(node)%values(1,-5:n(1)+6,-5:n(2)+6,-5:n(3)+6), &
        source_metrics%fields(node)%values(1,n(1)+6,n(2)+6,n(3)+6), &
        rates%fields(node)%values(nv,n(1),n(2),n(3)),stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      end if
      metrics%fields(node)%values(1,:,:,:)=extended_geometry(node)%volume_fraction
      source_metrics%fields(node)%values(1,:,:,:)=extended_geometry(node)%volume_fraction(4:n(1)+9,4:n(2)+9,4:n(3)+9)
      rates%fields(node)%values=0
    end do
    call agree(valid,ready)
    if(.not.ready) return
    call fill_mpi_eb_tree_halo_3d(species,tree,owners,geometry,state,temperature,ghost,gt,comm,ready, &
      inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
      boundary_faces=boundary_faces,inflow_state=inflow_state,ghost_volume_fraction=vf,halo_width=6, &
      target_volume_fraction=metrics,restrict_snapshot=.false., &
      ancestor_start=ancestor_start,ancestor_end=ancestor_end,temperature_start=temperature_start, &
      temperature_end=temperature_end,ancestor_fraction=ancestor_fraction)
    if(.not.ready) return
    faces='outflow'
    if(present(boundary_faces)) faces=boundary_faces
    if(present(inflow_mach)) then
      where(inflow_mach>=0) faces='inflow'
    end if
    call fill_sdc_mpi_eb_tree_source_halo_3d(tree,owners,geometry,predictor,nv,3,source,comm,ready, &
      boundary_faces=faces,target_volume_fraction=source_metrics, &
      ancestor_start=source_start,ancestor_end=source_end,ancestor_fraction=source_fraction)
    if(.not.ready) return
    if(present(predictor_halo)) then
      valid=allocated(predictor_halo%fields)
      if(valid) valid=size(predictor_halo%fields)==size(owners)
      call agree(valid,ready);if(.not.ready) return
      do node=1,size(owners)
        valid=allocated(predictor_halo%fields(node)%values).eqv.(owners(node)==rank)
        if(.not.valid) exit
        if(owners(node)/=rank) cycle
        valid=all(shape(predictor_halo%fields(node)%values)==shape(source%fields(node)%values))
        if(.not.valid) exit
        valid=all(ieee_is_finite(predictor_halo%fields(node)%values))
        if(.not.valid) exit
        n=tree%nodes(node)%shape
        offset=[lbound(predictor_halo%fields(node)%values,2),lbound(predictor_halo%fields(node)%values,3), &
          lbound(predictor_halo%fields(node)%values,4)]-1
        do k=1,n(3)
        do j=1,n(2)
        do i=1,n(1)
          if(geometry(node)%volume_fraction(i,j,k)<=0) cycle
          valid=valid.and.all(predictor_halo%fields(node)%values(:,i+3+offset(1),j+3+offset(2),k+3+offset(3))== &
            source%fields(node)%values(:,i+3,j+3,k+3))
        end do
        end do
        end do
        if(.not.valid) exit
        source%fields(node)%values=predictor_halo%fields(node)%values
      end do
      call agree(valid,ready);if(.not.ready) return
    end if
    do node=1,size(owners)
      if(owners(node)/=rank) cycle
      if(selected_depth>=0.and.tree%nodes(node)%level/=selected_depth) cycle
      n=tree%nodes(node)%shape
      call extract_eb_geometry_3d(extended_geometry(node),[4,4,4],n+9,working,valid)
      if(.not.valid) exit
      allocate(flat(0:n(1)+7,0:n(2)+7,0:n(3)+7),stat=ios)
      if(ios/=0) then
        valid=.false.;exit
      end if
      call compute_reactive_eb_flattening_3d(species,ghost%fields(node)%values,gt%fields(node)%values(1,:,:,:), &
        vf%fields(node)%values(1,:,:,:),n+6,flat,valid)
      if(.not.valid) exit
      physical_boundary%lower=4-int(origin(:,node))
      physical_boundary%upper=3+int(domain(:,tree%nodes(node)%level+1)-origin(:,node))
      physical_boundary%faces=faces
      work_stamp=tree_work_start()
      call compute_sdc_eb_hydro_3d(species,working,ghost%fields(node)%values,gt%fields(node)%values(1,:,:,:), &
        vf%fields(node)%values(1,:,:,:),dt,flat,source%fields(node)%values,rates%fields(node)%values, &
        candidate_flux(node),boundary,wall,valid,difmag=coefficient,boundary=physical_boundary)
      call tree_work_stop(node,work_stamp)
      if(.not.valid) exit
      deallocate(flat)
    end do
    call agree(valid,ready)
    if(.not.ready) return
    call move_alloc(rates%fields,hydro_rate%fields)
    call move_alloc(candidate_flux,flux);ok=.true.
  contains
    subroutine agree(condition,result)
      logical, intent(in) :: condition
      logical, intent(out) :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      if(ierr/=MPI_SUCCESS) result=.false.
    end subroutine
  end subroutine compute_sdc_mpi_eb_tree_hydro_3d
end module mpi_eb_tree_hydro_3d_mod
