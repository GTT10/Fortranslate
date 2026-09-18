module mpi_amr_patch_tree_transport_subcycle_3d_mod
  use mpi_f08
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use mpi_amr_patch_tree_chemistry_3d_mod, only: clone_mpi_amr_tree_inputs_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: recover_mpi_amr_tree_temperatures_3d
  use mpi_amr_patch_tree_reflux_3d_mod, only: reflux_mpi_amr_tree_3d
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use gas_transport_mod, only: gas_transport_species
  use reactive_transport_3d_mod, only: reactive_transport_timestep_3d,reactive_isothermal_timestep_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use mpi_amr_patch_tree_transport_level_3d_mod, only: advance_mpi_amr_tree_transport_level_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  implicit none
  private
  public :: compute_mpi_amr_tree_subcycle_transport_timestep_3d,advance_mpi_amr_tree_subcycle_transport_3d
contains
  subroutine compute_mpi_amr_tree_subcycle_transport_timestep_3d(tree,owners,species,transport,state,temperature, &
      root_spacing,cfl,viscosity,conduction,diffusion,dt,maximum_diffusivity,comm,ok,faces,wall_temperature)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: root_spacing(3),cfl
    type(gas_transport_species), intent(in) :: transport(:)
    logical, intent(in) :: viscosity,conduction,diffusion
    real(dp), intent(inout) :: dt,maximum_diffusivity
    logical, intent(out) :: ok
    real(dp), allocatable :: h(:,:)
    character(len=*), optional, intent(in) :: faces(6)
    real(dp), optional, intent(in) :: wall_temperature(:)
    real(dp) :: candidate,local_dt,local_diffusivity,diffusivity,wall_dt
    integer(int64), allocatable :: origin(:,:),domain(:,:),counts(:)
    integer, allocatable :: ratios(:)
    character(len=32) :: local_faces(6)
    integer :: node,n(3),axis,level,rank,status,ios
    integer :: flags(5),min_flags(5),max_flags(5)
    real(dp) :: controls(10),low(10),high(10),global_dt,global_diffusivity
    character(len=64) :: options(6),reference_options(6)
    logical :: valid
    ok=.false.
    call validate_mpi_amr_tree_local_3d(tree,owners,state,reactive_nvar(size(species)),comm,valid)
    if (.not.valid) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,valid)
    if (.not.valid) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    controls=0;controls(:4)=[root_spacing,cfl]
    valid=all(ieee_is_finite(controls)).and.all(controls(:4)>0).and.size(transport)==size(species)
    if (present(faces)) then
      if (.not.valid_reactive_boundaries_3d(faces)) valid=.false.
      if (any(faces=='isothermal_wall')) then
        if (.not.conduction.or..not.present(wall_temperature)) valid=.false.
      end if
    end if
    if (present(wall_temperature)) then
      if (size(wall_temperature)/=6.or..not.all(ieee_is_finite(wall_temperature))) then
        valid=.false.
      else
        controls(5:)=wall_temperature
      end if
    end if
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      if (.not.all(ieee_is_finite(state%fields(node)%values))) valid=.false.
      if (.not.all(ieee_is_finite(temperature%fields(node)%values))) valid=.false.
      if (any(temperature%fields(node)%values<=0)) valid=.false.
    end do
    call agree_dt(valid)
    if (.not.valid) return
    flags=[merge(1,0,viscosity),merge(1,0,conduction),merge(1,0,diffusion), &
      merge(1,0,present(faces)),merge(1,0,present(wall_temperature))]
    call MPI_Allreduce(flags,min_flags,5,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,max_flags,5,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(min_flags/=max_flags)) return
    call MPI_Allreduce(controls,low,10,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,10,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(low/=high)) return
    options='outflow'
    if (present(faces)) options=faces
    reference_options=options
    call MPI_Bcast(reference_options,64*6,MPI_CHARACTER,0,comm,status)
    if (status/=MPI_SUCCESS) return
    valid=all(options==reference_options);call agree_dt(valid)
    if (.not.valid) return
    allocate(h(3,size(tree%nodes)),stat=ios)
    valid=ios==0;call agree_dt(valid)
    if (.not.valid) return
    call layout(tree,owners,species,state,temperature,root_spacing,origin,domain,ratios,counts,comm,valid)
    if (.not.valid) return
    do node=1,size(tree%nodes)
      h(:,node)=root_spacing*real(domain(:,1),dp)/real(domain(:,tree%nodes(node)%level+1),dp)
    end do
    candidate=huge(1.0_dp);diffusivity=0
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape;level=tree%nodes(node)%level+1
      call reactive_transport_timestep_3d(species,transport,state%fields(node)%values,temperature%fields(node)%values(1,:,:,:), &
        n(1),n(2),n(3),h(1,node),h(2,node),h(3,node),cfl,viscosity,conduction,diffusion,local_dt,local_diffusivity,valid)
      if (.not.valid) exit
      if (local_dt<=huge(1._dp)/real(counts(level),dp)) candidate=min(candidate,local_dt*real(counts(level),dp))
      diffusivity=max(diffusivity,local_diffusivity)
      if (present(faces)) then
        if (any(faces=='isothermal_wall')) then
          local_faces='outflow';level=tree%nodes(node)%level+1
          do axis=1,3
            if (origin(axis,node)==0.and.faces(2*axis-1)=='isothermal_wall') &
              local_faces(2*axis-1)='isothermal_wall'
            if (origin(axis,node)+int(n(axis),int64)==domain(axis,level).and.faces(2*axis)=='isothermal_wall') &
              local_faces(2*axis)='isothermal_wall'
          end do
          call reactive_isothermal_timestep_3d(species,transport,state%fields(node)%values, &
            temperature%fields(node)%values(1,:,:,:),h(:,node),local_faces,wall_temperature,cfl,wall_dt,valid)
          if (.not.valid) exit
          if (wall_dt<=huge(1._dp)/real(counts(level),dp)) candidate=min(candidate,wall_dt*real(counts(level),dp))
        end if
      end if
    end do
    call agree_dt(valid)
    if (.not.valid) return
    call MPI_Allreduce(candidate,global_dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(diffusivity,global_diffusivity,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    dt=global_dt;maximum_diffusivity=global_diffusivity;ok=.true.
  contains
    subroutine agree_dt(condition)
      logical, intent(inout) :: condition
      logical :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      condition=.false.
      if (ierr==MPI_SUCCESS) condition=result
    end subroutine
  end subroutine
  subroutine advance_mpi_amr_tree_subcycle_transport_3d(tree,owners,species,transport,state,temperature,spacing,dt,faces, &
      viscosity,conduction,diffusion,barodiffusion,boundary_change,maximum_correction,minimum_theta,comm,ok, &
      inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter,level_steps,outlet_pressure)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: spacing(3),dt
    type(gas_transport_species), intent(in) :: transport(:)
    character(len=*), intent(in) :: faces(6)
    character(len=*), optional, intent(in) :: coarse_limiter
    logical, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), optional, intent(in) :: wall_temperature(:),wall_velocity(:,:)
    real(dp), intent(inout) :: boundary_change(:),maximum_correction,minimum_theta
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:)
    real(dp), optional, intent(in) :: outlet_pressure(6)
    integer(int64), optional, intent(inout) :: level_steps(:)
    type(amr_tree_state_3d) :: candidate,temps,old_endpoint,new_endpoint
    type(amr_tree_flux_3d), allocatable :: integrated(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:),counts(:),executed(:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: starts(:),finishes(:)
    real(dp) :: budget(size(boundary_change)),largest,theta,dt_low,dt_high
    integer :: nvar,node,n(3),rank,status,ios
    logical :: valid
    ok=.false.
    call layout(tree,owners,species,state,temperature,spacing,origin,domain,ratios,counts,comm,valid)
    if (.not.valid) return
    nvar=reactive_nvar(size(species))
    valid=size(boundary_change)==nvar.and.ieee_is_finite(dt)
    if (valid) valid=dt>0
    if (present(level_steps)) valid=valid.and.size(level_steps)==size(counts)
    call agree(comm,valid)
    if (.not.valid) return
    call MPI_Allreduce(dt,dt_low,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(dt,dt_high,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (dt_low/=dt_high) return
    call clone_mpi_amr_tree_inputs_3d(tree,owners,size(species),state,temperature, &
      candidate,temps,comm,rank,valid)
    if (.not.valid) return
    call average_down_mpi_amr_tree_3d(tree,owners,candidate,nvar,comm,valid)
    if (.not.valid) return
    call recover_from(0,valid)
    if (.not.valid) return
    allocate(old_endpoint%fields(size(tree%nodes)),new_endpoint%fields(size(tree%nodes)),stat=ios)
    valid=ios==0;call agree(comm,valid)
    if (.not.valid) return
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      allocate(old_endpoint%fields(node)%values(nvar+1,n(1),n(2),n(3)),stat=ios)
      if (ios/=0) exit
      allocate(new_endpoint%fields(node)%values(nvar+1,n(1),n(2),n(3)),stat=ios)
      if (ios/=0) exit
    end do
    valid=ios==0;call agree(comm,valid)
    if (.not.valid) return
    allocate(starts(size(counts)),finishes(size(counts)),executed(size(counts)),stat=ios)
    valid=ios==0;call agree(comm,valid)
    if (.not.valid) return
    starts=0;finishes=1;executed=0;budget=0;largest=0;theta=1
    call interval(0,0._dp,1._dp,integrated,valid)
    if (.not.valid) return
    if (any(executed/=counts)) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(temps%fields,temperature%fields)
    boundary_change=budget;maximum_correction=largest;minimum_theta=theta
    if (present(level_steps)) level_steps=executed
    ok=.true.
  contains
    recursive subroutine interval(depth,begin_time,end_time,effective,success)
      integer, intent(in) :: depth
      real(dp), intent(in) :: begin_time,end_time
      type(amr_tree_flux_3d), allocatable, intent(out) :: effective(:)
      logical, intent(out) :: success
      type(amr_tree_flux_3d), allocatable :: total(:),child_flux(:)
      real(dp) :: fraction0(depth),fraction1(depth),correction(nvar),max_change,sub_begin,sub_end,local_dt,local_theta
      integer :: peer,ancestor,allocation_status
      integer(int64) :: substep,ratio
      logical :: good
      success=.false.
      if (end_time<=begin_time) return
      local_dt=dt*(end_time-begin_time)
      if (local_dt<=0) return
      starts(depth+1)=begin_time;finishes(depth+1)=end_time
      do ancestor=1,depth
        fraction0(ancestor)=max(0._dp,min(1._dp,(begin_time-starts(ancestor))/(finishes(ancestor)-starts(ancestor))))
        fraction1(ancestor)=max(0._dp,min(1._dp,(end_time-starts(ancestor))/(finishes(ancestor)-starts(ancestor))))
      end do
      call pack_depth(depth,old_endpoint)
      local_theta=1
      call advance_mpi_amr_tree_transport_level_3d(tree,owners,species,transport,candidate,temps,depth,spacing,local_dt, &
        faces,viscosity,conduction,diffusion,barodiffusion,total,local_theta,comm,good, &
        old_endpoint,new_endpoint,fraction0,fraction1,inflow_state,inflow_temperature, &
        wall_temperature,wall_velocity,coarse_limiter,outlet_pressure)
      if (.not.good) return
      theta=min(theta,local_theta)
      executed(depth+1)=executed(depth+1)+1
      call pack_depth(depth,new_endpoint)
      if (depth==0) then
        if (rank==owners(1)) call root_budget(total(1))
        call MPI_Bcast(budget,nvar,MPI_DOUBLE_PRECISION,owners(1),comm,status)
        if (status/=MPI_SUCCESS) return
      end if
      if (depth+1<size(counts)) then
        ratio=int(ratios(depth+2),int64)**2
        sub_begin=begin_time
        do substep=1,ratio
          sub_end=begin_time+(end_time-begin_time)*real(substep,dp)/real(ratio,dp)
          if (substep==ratio) sub_end=end_time
          call interval(depth+1,sub_begin,sub_end,child_flux,good)
          if (.not.good) return
          allocation_status=0
          do peer=1,size(tree%nodes)
            if (owners(peer)/=rank) cycle
            if (tree%nodes(peer)%level<=depth) cycle
            if (.not.allocated(total(peer)%x)) then
              allocate(total(peer)%x,source=child_flux(peer)%x,stat=allocation_status)
              if (allocation_status/=0) exit
              allocate(total(peer)%y,source=child_flux(peer)%y,stat=allocation_status)
              if (allocation_status/=0) exit
              allocate(total(peer)%z,source=child_flux(peer)%z,stat=allocation_status)
              if (allocation_status/=0) exit
            else
              total(peer)%x=total(peer)%x+child_flux(peer)%x
              total(peer)%y=total(peer)%y+child_flux(peer)%y
              total(peer)%z=total(peer)%z+child_flux(peer)%z
            end if
          end do
          good=allocation_status==0;call agree(comm,good)
          if (.not.good) return
          sub_begin=sub_end
        end do
      end if
      correction=0;max_change=0
      call reflux_mpi_amr_tree_3d(tree,owners,candidate,nvar,total,spacing,faces,correction,max_change,comm,good, &
        minimum_depth=depth,synchronized_flux=effective)
      if (.not.good) return
      largest=max(largest,max_change)
      if (depth==0) budget=budget+correction
      call recover_from(depth,good)
      if (.not.good) return
      success=.true.
    end subroutine
    subroutine pack_depth(depth,packed)
      integer, intent(in) :: depth
      type(amr_tree_state_3d), intent(inout) :: packed
      integer :: peer
      do peer=1,size(tree%nodes)
        if (owners(peer)/=rank.or.tree%nodes(peer)%level/=depth) cycle
        packed%fields(peer)%values(:nvar,:,:,:)=candidate%fields(peer)%values
        packed%fields(peer)%values(nvar+1,:,:,:)=temps%fields(peer)%values(1,:,:,:)
      end do
    end subroutine
    subroutine recover_from(depth,success)
      integer, intent(in) :: depth
      logical, intent(out) :: success
      integer :: level
      do level=depth,tree%level_count()-1
        call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,candidate,temps,comm,success,depth=level)
        if (.not.success) return
      end do
      success=.true.
    end subroutine
    subroutine root_budget(flux)
      type(amr_tree_flux_3d), intent(in) :: flux
      integer :: i,j,k,dims(3)
      dims=tree%nodes(1)%shape
      if (faces(1)/='periodic') then
        do k=1,dims(3)
          do j=1,dims(2)
            budget=budget+(flux%x(:,1,j,k)-flux%x(:,dims(1)+1,j,k))*spacing(2)*spacing(3)
          end do
        end do
      end if
      if (faces(3)/='periodic') then
        do k=1,dims(3)
          do i=1,dims(1)
            budget=budget+(flux%y(:,i,1,k)-flux%y(:,i,dims(2)+1,k))*spacing(1)*spacing(3)
          end do
        end do
      end if
      if (faces(5)/='periodic') then
        do j=1,dims(2)
          do i=1,dims(1)
            budget=budget+(flux%z(:,i,j,1)-flux%z(:,i,j,dims(3)+1))*spacing(1)*spacing(2)
          end do
        end do
      end if
    end subroutine
  end subroutine

  subroutine layout(tree,owners,species,state,temperature,spacing,origin,domain,ratios,counts,comm,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: spacing(3)
    integer(int64), allocatable, intent(out) :: origin(:,:),domain(:,:),counts(:)
    integer, allocatable, intent(out) :: ratios(:)
    logical, intent(out) :: ok
    integer :: level,node,rank,status,ios
    integer(int64) :: factor
    real(dp) :: low(3),high(3)
    logical :: valid
    call validate_mpi_amr_tree_local_3d(tree,owners,state,reactive_nvar(size(species)),comm,ok)
    if (.not.ok) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,ok)
    if (.not.ok) return
    ok=.false.
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    valid=all(ieee_is_finite(spacing))
    if (valid) valid=all(spacing>0)
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      if (.not.all(ieee_is_finite(state%fields(node)%values))) valid=.false.
      if (.not.all(ieee_is_finite(temperature%fields(node)%values))) then
        valid=.false.
      else if (any(temperature%fields(node)%values<=0)) then
        valid=.false.
      end if
    end do
    call agree(comm,valid)
    if (.not.valid) return
    call MPI_Allreduce(spacing,low,3,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(spacing,high,3,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(low/=high)) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
    call agree(comm,ok)
    if (.not.ok) return
    ok=.false.
    allocate(counts(size(ratios)),stat=ios)
    valid=ios==0;call agree(comm,valid)
    if (.not.valid) return
    counts(1)=1
    do level=2,size(ratios)
      factor=int(ratios(level),int64)**2
      if (counts(level-1)>huge(1_int64)/factor) return
      counts(level)=counts(level-1)*factor
    end do
    ok=.true.
  end subroutine
  subroutine agree(comm,condition)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(inout) :: condition
    logical :: result
    integer :: status
    call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,status)
    condition=.false.
    if (status==MPI_SUCCESS) condition=result
  end subroutine
end module mpi_amr_patch_tree_transport_subcycle_3d_mod
