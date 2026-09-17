module mpi_amr_patch_tree_transport_3d_mod
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use gas_transport_mod, only: gas_transport_species
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use reactive_transport_3d_mod, only: physical_transport_ghosted_fluxes_3d, &
    reactive_transport_timestep_3d,reactive_isothermal_timestep_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use mpi_amr_patch_tree_ghost_3d_mod, only: fill_mpi_amr_tree_ghosts_3d
  use mpi_amr_patch_tree_reflux_3d_mod, only: reflux_mpi_amr_tree_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: recover_mpi_amr_tree_temperatures_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  implicit none
  private
  public :: advance_mpi_amr_tree_transport_3d,compute_mpi_amr_tree_transport_timestep_3d
contains
  subroutine compute_mpi_amr_tree_transport_timestep_3d(tree,owners,species,transport,state,temperature, &
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
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    character(len=32) :: local_faces(6)
    integer :: node,n(3),axis,level,rank,status,ios,parent
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
    h(:,1)=root_spacing
    do node=2,size(tree%nodes)
      parent=tree%nodes(node)%parent
      h(:,node)=h(:,parent)/real(tree%nodes(node)%patch%refinement_ratio,dp)
      if (any(h(:,node)<=0)) return
    end do
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree_dt(valid)
    if (.not.valid) return
    candidate=huge(1.0_dp);diffusivity=0
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      call reactive_transport_timestep_3d(species,transport,state%fields(node)%values,temperature%fields(node)%values(1,:,:,:), &
        n(1),n(2),n(3),h(1,node),h(2,node),h(3,node),cfl,viscosity,conduction,diffusion,local_dt,local_diffusivity,valid)
      if (.not.valid) exit
      candidate=min(candidate,local_dt);diffusivity=max(diffusivity,local_diffusivity)
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
          candidate=min(candidate,wall_dt)
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
  subroutine advance_mpi_amr_tree_transport_3d(tree,owners,species,transport,state,temperature,root_spacing,dt,faces, &
      viscosity,conduction,diffusion,barodiffusion,boundary_change,maximum_correction,minimum_theta,comm,ok, &
      inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: root_spacing(3),dt
    type(gas_transport_species), intent(in) :: transport(:)
    character(len=*), intent(in) :: faces(6)
    logical, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(inout) :: boundary_change(:),maximum_correction,minimum_theta
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    character(len=*), optional, intent(in) :: coarse_limiter
    type(amr_tree_state_3d) :: initial,candidate,temps,initial_temperature
    real(dp), allocatable :: h(:,:),inlet(:,:)
    real(dp) :: first_budget(size(boundary_change)),second_budget(size(boundary_change)), &
      first_max,second_max,first_theta,second_theta
    integer :: node,nvar,rank,status,ios,parent
    real(dp) :: controls(4),low(4),high(4)
    real(dp), allocatable :: reference_inlet(:,:)
    character(len=64) :: options(7),reference_options(7)
    integer :: flags(7),min_flags(7),max_flags(7)
    real(dp) :: walls(24),reference_walls(24)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    logical :: valid
    ok=.false.
    nvar=reactive_nvar(size(species))
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nvar,comm,valid)
    if (.not.valid) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,valid)
    if (.not.valid) return
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    controls=[root_spacing,dt]
    valid=all(ieee_is_finite(controls)).and.all(controls>0)
    valid=valid.and.size(boundary_change)==nvar.and.valid_reactive_boundaries_3d(faces)
    valid=valid.and..not.(barodiffusion.and..not.diffusion).and.size(transport)==size(species)
    if (present(coarse_limiter)) valid=valid.and.(coarse_limiter=='mc'.or.coarse_limiter=='minmod')
    if (present(wall_temperature)) then
      if (size(wall_temperature)/=6.or..not.all(ieee_is_finite(wall_temperature))) valid=.false.
    end if
    if (present(wall_velocity)) then
      if (any(shape(wall_velocity)/=[3,6]).or..not.all(ieee_is_finite(wall_velocity))) valid=.false.
    end if
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      if (.not.all(ieee_is_finite(temperature%fields(node)%values))) valid=.false.
      if (any(temperature%fields(node)%values<=0)) valid=.false.
    end do
    call agree(valid)
    if (.not.valid) return
    call MPI_Allreduce(controls,low,4,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,4,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(low/=high)) return
    flags=[merge(1,0,viscosity),merge(1,0,conduction),merge(1,0,diffusion),merge(1,0,barodiffusion), &
      merge(1,0,present(wall_temperature)),merge(1,0,present(wall_velocity)),merge(1,0,present(coarse_limiter))]
    call MPI_Allreduce(flags,min_flags,7,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,max_flags,7,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(min_flags/=max_flags)) return
    options(:6)=faces;options(7)='pcm'
    if (present(coarse_limiter)) options(7)=coarse_limiter
    reference_options=options
    call MPI_Bcast(reference_options,64*7,MPI_CHARACTER,0,comm,status)
    if (status/=MPI_SUCCESS) return
    valid=all(options==reference_options)
    call agree(valid)
    if (.not.valid) return
    if (any(faces=='inflow')) then
      valid=present(inflow_state).and.present(inflow_temperature)
      if (valid) then
        valid=all(shape(inflow_state)==[nvar,6]).and.size(inflow_temperature)==6
        if (valid) valid=all(ieee_is_finite(inflow_state)).and.all(ieee_is_finite(inflow_temperature))
      end if
      call agree(valid)
      if (.not.valid) return
      allocate(inlet(nvar+1,6),reference_inlet(nvar+1,6),stat=ios)
      valid=ios==0;call agree(valid)
      if (.not.valid) return
      inlet(:nvar,:)=inflow_state;inlet(nvar+1,:)=inflow_temperature
      reference_inlet=inlet
      call MPI_Bcast(reference_inlet,size(reference_inlet),MPI_DOUBLE_PRECISION,0,comm,status)
      if (status/=MPI_SUCCESS) return
      valid=all(inlet==reference_inlet);call agree(valid)
      if (.not.valid) return
    end if
    walls=0
    if (present(wall_temperature)) walls(:6)=wall_temperature
    if (present(wall_velocity)) walls(7:)=reshape(wall_velocity,[18])
    reference_walls=walls
    call MPI_Bcast(reference_walls,24,MPI_DOUBLE_PRECISION,0,comm,status)
    if (status/=MPI_SUCCESS) return
    valid=all(walls==reference_walls);call agree(valid)
    if (.not.valid) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    call agree(valid)
    if (.not.valid) return
    allocate(h(3,size(tree%nodes)),stat=ios)
    valid=ios==0;call agree(valid)
    if (.not.valid) return
    h(:,1)=root_spacing
    do node=2,size(tree%nodes)
      parent=tree%nodes(node)%parent
      h(:,node)=h(:,parent)/real(tree%nodes(node)%patch%refinement_ratio,dp)
      if (any(h(:,node)<=0)) return
    end do
    call copy_local(state,initial,valid)
    if (.not.valid) return
    call copy_local(temperature,temps,valid)
    if (.not.valid) return
    call average_down_mpi_amr_tree_3d(tree,owners,initial,nvar,comm,valid)
    if (.not.valid) return
    call recover_tree(initial,temps,valid)
    if (.not.valid) return
    call copy_local(initial,candidate,valid)
    if (.not.valid) return
    call copy_local(temps,initial_temperature,valid)
    if (.not.valid) return
    call euler_stage(candidate,temps,first_budget,first_max,first_theta,valid)
    if (.not.valid) return
    call euler_stage(candidate,temps,second_budget,second_max,second_theta,valid)
    if (.not.valid) return
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      candidate%fields(node)%values=.5_dp*initial%fields(node)%values+.5_dp*candidate%fields(node)%values
    end do
    call average_down_mpi_amr_tree_3d(tree,owners,candidate,nvar,comm,valid)
    if (.not.valid) return
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      temps%fields(node)%values=.5_dp*(initial_temperature%fields(node)%values+temps%fields(node)%values)
    end do
    call recover_tree(candidate,temps,valid)
    if (.not.valid) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(temps%fields,temperature%fields)
    boundary_change=.5_dp*(first_budget+second_budget)
    maximum_correction=max(first_max,second_max);minimum_theta=min(first_theta,second_theta);ok=.true.
  contains
    subroutine recover_tree(u,t,success)
      type(amr_tree_state_3d), intent(in) :: u
      type(amr_tree_state_3d), intent(inout) :: t
      logical, intent(out) :: success
      call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,u,t,comm,success)
    end subroutine
    subroutine euler_stage(u,t,budget,largest,theta,success)
      type(amr_tree_state_3d), intent(inout) :: u,t
      real(dp), intent(out) :: budget(:),largest,theta
      logical, intent(out) :: success
      type(amr_tree_state_3d) :: packed,updated,ghosts
      type(amr_tree_flux_3d), allocatable :: flux(:)
      real(dp), allocatable :: fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:)
      real(dp) :: correction(nvar),rhs(nvar),local_theta,global_theta
      character(len=32) :: local_faces(6)
      integer :: peer,n(3),dims(3),i,j,k,allocation_status,axis,level
      success=.false.;budget=0;largest=0;correction=0;theta=1
      allocate(packed%fields(size(tree%nodes)),flux(size(tree%nodes)),stat=allocation_status)
      success=allocation_status==0;call agree(success)
      if (.not.success) return
      do peer=1,size(tree%nodes)
        if (owners(peer)/=rank) cycle
        n=tree%nodes(peer)%shape
        allocate(packed%fields(peer)%values(nvar+1,n(1),n(2),n(3)),stat=allocation_status)
        if (allocation_status/=0) exit
        packed%fields(peer)%values(:nvar,:,:,:)=u%fields(peer)%values
        packed%fields(peer)%values(nvar+1,:,:,:)=t%fields(peer)%values(1,:,:,:)
      end do
      success=allocation_status==0;call agree(success)
      if (.not.success) return
      call copy_local(u,updated,success)
      if (.not.success) return
      call fill_mpi_amr_tree_ghosts_3d(tree,owners,packed,nvar+1,faces,ghosts,comm,success, &
        viscous_wall=.true.,inflow_state=inlet,wall_velocity=wall_velocity, &
        coarse_limiter=coarse_limiter,species_count=size(species))
      if (.not.success) return
      do peer=1,size(tree%nodes)
        if (owners(peer)/=rank) cycle
        n=tree%nodes(peer)%shape;dims=n+4
        associate(ghost=>ghosts%fields(peer)%values)
          allocate(fx(nvar,n(1)+1,n(2),n(3)),fy(nvar,n(1),n(2)+1,n(3)), &
            fz(nvar,n(1),n(2),n(3)+1),stat=allocation_status)
          if (allocation_status/=0) then
            success=.false.;exit
          end if
          local_faces='outflow';level=tree%nodes(peer)%level+1
          do axis=1,3
            if (origin(axis,peer)==0) local_faces(2*axis-1)=faces(2*axis-1)
            if (origin(axis,peer)+int(n(axis),int64)==domain(axis,level)) local_faces(2*axis)=faces(2*axis)
          end do
          call physical_transport_ghosted_fluxes_3d(species,transport,u%fields(peer)%values, &
            t%fields(peer)%values(1,:,:,:),ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:),n,h(:,peer),dt,local_faces, &
            viscosity,conduction,diffusion,barodiffusion,fx,fy,fz,local_theta,success,wall_temperature,wall_velocity)
          theta=min(theta,local_theta)
          if (.not.success) exit
          do k=1,n(3)
            do j=1,n(2)
              do i=1,n(1)
                rhs=-(fx(:,i+1,j,k)-fx(:,i,j,k))/h(1,peer) &
                  -(fy(:,i,j+1,k)-fy(:,i,j,k))/h(2,peer) &
                  -(fz(:,i,j,k+1)-fz(:,i,j,k))/h(3,peer)
                updated%fields(peer)%values(:,i,j,k)=u%fields(peer)%values(:,i,j,k)+dt*rhs
              end do
            end do
          end do
          allocate(flux(peer)%x,source=dt*fx,stat=allocation_status)
          if (allocation_status/=0) then
            success=.false.;exit
          end if
          allocate(flux(peer)%y,source=dt*fy,stat=allocation_status)
          if (allocation_status/=0) then
            success=.false.;exit
          end if
          allocate(flux(peer)%z,source=dt*fz,stat=allocation_status)
          if (allocation_status/=0) then
            success=.false.;exit
          end if
          deallocate(fx,fy,fz)
        end associate
      end do
      call agree(success)
      if (.not.success) return
      call MPI_Allreduce(theta,global_theta,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if (status/=MPI_SUCCESS) then
        success=.false.;return
      end if
      theta=global_theta
      if (rank==owners(1)) then
        n=tree%nodes(1)%shape
        do k=1,n(3)
          do j=1,n(2)
            budget=budget+(flux(1)%x(:,1,j,k)-flux(1)%x(:,n(1)+1,j,k))*h(2,1)*h(3,1)
          end do
        end do
        do k=1,n(3)
          do i=1,n(1)
            budget=budget+(flux(1)%y(:,i,1,k)-flux(1)%y(:,i,n(2)+1,k))*h(1,1)*h(3,1)
          end do
        end do
        do j=1,n(2)
          do i=1,n(1)
            budget=budget+(flux(1)%z(:,i,j,1)-flux(1)%z(:,i,j,n(3)+1))*h(1,1)*h(2,1)
          end do
        end do
      end if
      call MPI_Bcast(budget,nvar,MPI_DOUBLE_PRECISION,owners(1),comm,status)
      if (status/=MPI_SUCCESS) then
        success=.false.;return
      end if
      call reflux_mpi_amr_tree_3d(tree,owners,updated,nvar,flux,root_spacing,faces,correction,largest,comm,success)
      if (.not.success) return
      budget=budget+correction
      call recover_tree(updated,t,success)
      if (.not.success) return
      call move_alloc(updated%fields,u%fields)
    end subroutine
    subroutine agree(condition)
      logical, intent(inout) :: condition
      logical :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      condition=.false.
      if (ierr==MPI_SUCCESS) condition=result
    end subroutine
    subroutine copy_local(source,destination,success)
      type(amr_tree_state_3d), intent(in) :: source
      type(amr_tree_state_3d), intent(inout) :: destination
      logical, intent(out) :: success
      type(amr_tree_state_3d) :: copy
      integer :: peer,ierr
      allocate(copy%fields(size(tree%nodes)),stat=ierr)
      if (ierr==0) then
        do peer=1,size(tree%nodes)
          if (owners(peer)/=rank) cycle
          allocate(copy%fields(peer)%values,source=source%fields(peer)%values,stat=ierr)
          if (ierr/=0) exit
        end do
      end if
      success=ierr==0;call agree(success)
      if (success) call move_alloc(copy%fields,destination%fields)
    end subroutine
  end subroutine advance_mpi_amr_tree_transport_3d
end module
