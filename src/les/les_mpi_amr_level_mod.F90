module les_mpi_amr_level_mod
  use mpi_amr_tree_work_mod, only: tree_work_start,tree_work_stop
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use les_smagorinsky_3d_mod, only: compute_les_ghosted_rhs_3d
  use les_dynamic_smagorinsky_mod, only: compute_dynamic_les_rhs_ghosted_3d
  use les_filter_3d_mod, only: les_filter_weights
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use amr_patch_tree_ghost_3d_mod, only: apply_amr_tree_pressure_ghosts_3d
  use mpi_amr_patch_tree_ghost_3d_mod, only: fill_mpi_amr_tree_ghosts_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: recover_mpi_amr_tree_temperatures_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  implicit none
  private
  public :: advance_mpi_amr_tree_les_level_3d
contains
  subroutine advance_mpi_amr_tree_les_level_3d(tree,owners,species,state,temperature,depth,root_spacing,dt,faces, &
      cs,prt,sct,integrated_flux,maximum_diffusivity,comm,ok, &
      coarse_start,coarse_end,fractions_start,fractions_end,inflow_state,inflow_temperature, &
      wall_temperature,wall_velocity,coarse_limiter,outlet_pressure, &
      model,ci,test_model,test_ratio,grid_ratio,harmonic_mean,inflow_mach,inflow_pressure)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: root_spacing(3),dt
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: cs,prt,sct
    real(dp), optional, intent(in) :: ci
    character(len=*), optional, intent(in) :: model,test_model
    integer, optional, intent(in) :: test_ratio,grid_ratio
    logical, optional, intent(in) :: harmonic_mean
    real(dp) :: isotropic
    character(len=64) :: selected_model,filter
    integer :: filter_ratio,grid,width
    logical :: dynamic,harmonic,filter_valid
    real(dp), allocatable :: weights(:)
    real(dp), intent(inout) :: maximum_diffusivity
    integer, intent(in) :: depth
    type(amr_tree_flux_3d), allocatable, intent(inout) :: integrated_flux(:)
    type(amr_tree_flux_3d), allocatable :: first_flux(:),second_flux(:)
    type(amr_tree_state_3d), optional, intent(in) :: coarse_start,coarse_end
    real(dp), optional, intent(in) :: fractions_start(:),fractions_end(:)
    integer :: endpoint_flags(5),minimum_flags(5),maximum_flags(5)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6)
    real(dp) :: pressure(6),mach(6),tin(6),pin(6)
    character(len=*), optional, intent(in) :: coarse_limiter
    type(amr_tree_state_3d) :: initial,candidate,temps,first_temperature
    real(dp), allocatable :: h(:,:),inlet(:,:)
    real(dp) :: first_diffusivity,second_diffusivity
    integer :: node,nvar,rank,status,ios,parent
    real(dp) :: controls(35),low(35),high(35)
    real(dp), allocatable :: reference_inlet(:,:)
    character(len=64) :: options(9),reference_options(9)
    integer :: flags(3),min_flags(3),max_flags(3)
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
    pressure=0
    if (present(outlet_pressure)) pressure=outlet_pressure
    isotropic=0
    if (present(ci)) isotropic=ci
    selected_model='smagorinsky'
    if (present(model)) selected_model=model
    dynamic=selected_model=='dynamic';filter='box';filter_ratio=2;grid=1;harmonic=.false.;width=2
    if(present(test_model)) filter=test_model
    if(present(test_ratio)) filter_ratio=test_ratio
    if(present(grid_ratio)) grid=grid_ratio
    if(present(harmonic_mean)) harmonic=harmonic_mean
    filter_valid=.true.
    if(dynamic) then
      call les_filter_weights(filter,filter_ratio,weights,filter_valid)
      if(filter_valid) width=(size(weights)-1)/2+5
      filter_valid=filter_valid.and.grid>=1.and.grid<=64
      filter_valid=filter_valid.and.all(faces=='periodic'.or.faces=='inflow'.or.faces=='outflow'.or. &
        faces=='slip_wall'.or.faces=='no_slip_wall'.or.faces=='isothermal_wall')
    else
      if(present(test_model).or.present(test_ratio).or.present(grid_ratio).or.present(harmonic_mean)) filter_valid=.false.
    end if
    mach=-1;pin=0;tin=0
    if(present(inflow_mach)) mach=inflow_mach
    if(present(inflow_pressure)) pin=inflow_pressure
    valid=.true.
    if(present(inflow_temperature)) then
      valid=size(inflow_temperature)==6
      if(valid) then
        where(mach>=0) tin=inflow_temperature
      end if
    end if
    valid=valid.and.(.not.present(inflow_pressure).or.present(inflow_mach))
    valid=valid.and.all(ieee_is_finite(mach)).and.all(mach<=100)
    valid=valid.and..not.any(mach>=0.and.faces/='inflow')
    valid=valid.and.all(ieee_is_finite(pin)).and.all(pin>=0).and.all(pin<=1e15_dp)
    valid=valid.and..not.any(pin>0.and.mach<0)
    valid=valid.and.all(ieee_is_finite(tin)).and..not.any(tin<=0.and.mach>=0)
    call agree(valid)
    if(.not.valid) return
    controls=[root_spacing,dt,pressure,cs,prt,sct,isotropic,real(filter_ratio,dp),real(grid,dp), &
      merge(1._dp,0._dp,harmonic),mach,tin,pin]
    valid=all(ieee_is_finite(controls)).and.all(controls(:4)>0)
    valid=valid.and.all(pressure>=0).and.all(pressure<=1e15_dp)
    do node=1,6
      if (pressure(node)>0.and.faces(node)/='outflow') valid=.false.
    end do
    valid=valid.and.valid_reactive_boundaries_3d(faces)
    valid=valid.and.dt<=1e6_dp.and.cs>=0.and.cs<=2.and.min(prt,sct)>=1e-6_dp.and.max(prt,sct)<=1e6_dp
    valid=valid.and.isotropic>=0.and.isotropic<=1.and.filter_valid
    valid=valid.and.(selected_model=='smagorinsky'.or.selected_model=='wale'.or.selected_model=='vreman'.or.dynamic)
    if (present(coarse_limiter)) valid=valid.and.(coarse_limiter=='mc'.or.coarse_limiter=='minmod')
    if (present(wall_temperature)) then
      if (size(wall_temperature)/=6.or..not.all(ieee_is_finite(wall_temperature))) valid=.false.
    end if
    if (present(wall_velocity)) then
      if (any(shape(wall_velocity)/=[3,6]).or..not.all(ieee_is_finite(wall_velocity))) valid=.false.
    end if
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      if (.not.all(ieee_is_finite(state%fields(node)%values))) valid=.false.
      if (.not.all(ieee_is_finite(temperature%fields(node)%values))) valid=.false.
      if (any(temperature%fields(node)%values<=0)) valid=.false.
    end do
    call agree(valid)
    if (.not.valid) return
    call MPI_Allreduce(controls,low,35,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,35,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(low/=high)) return
    flags=[merge(1,0,present(wall_temperature)),merge(1,0,present(wall_velocity)),merge(1,0,present(coarse_limiter))]
    call MPI_Allreduce(flags,min_flags,3,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,max_flags,3,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(min_flags/=max_flags)) return
    options(:6)=faces;options(7)='pcm';options(8)=selected_model;options(9)=filter
    if (present(coarse_limiter)) options(7)=coarse_limiter
    reference_options=options
    call MPI_Bcast(reference_options,64*9,MPI_CHARACTER,0,comm,status)
    if (status/=MPI_SUCCESS) return
    valid=all(options==reference_options)
    call agree(valid)
    if (.not.valid) return
    endpoint_flags=[depth,merge(1,0,present(coarse_start)),merge(1,0,present(coarse_end)), &
      merge(1,0,present(fractions_start)),merge(1,0,present(fractions_end))]
    call MPI_Allreduce(endpoint_flags,minimum_flags,5,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(endpoint_flags,maximum_flags,5,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(minimum_flags/=maximum_flags)) return
    if (depth<0.or.depth>=tree%level_count().or.any(endpoint_flags(2:)/=endpoint_flags(2))) return
    if (depth>0.and.endpoint_flags(2)==0) return
    if (present(coarse_start)) then
      valid=size(fractions_start)==depth.and.size(fractions_end)==depth
      if (valid) valid=all(ieee_is_finite(fractions_start)).and.all(ieee_is_finite(fractions_end))
      if (valid) valid=all(fractions_start>=0).and.all(fractions_end<=1).and.all(fractions_end>=fractions_start)
      call agree(valid)
      if (.not.valid) return
      call validate_mpi_amr_tree_local_3d(tree,owners,coarse_start,nvar+1,comm,valid,maximum_depth=depth-1)
      if (.not.valid) return
      call validate_mpi_amr_tree_local_3d(tree,owners,coarse_end,nvar+1,comm,valid,maximum_depth=depth-1)
      if (.not.valid) return
      do node=1,size(tree%nodes)
        if (owners(node)/=rank.or.tree%nodes(node)%level>=depth) cycle
        if (.not.all(ieee_is_finite(coarse_start%fields(node)%values))) valid=.false.
        if (.not.all(ieee_is_finite(coarse_end%fields(node)%values))) valid=.false.
      end do
      call agree(valid)
      if (.not.valid) return
    end if
    if (any(faces=='inflow'.and.mach<0)) then
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
    call copy_local(initial,candidate,valid)
    if (.not.valid) return
    call euler_stage(candidate,temps,first_flux,1,first_diffusivity,valid)
    if (.not.valid) return
    call copy_local(temps,first_temperature,valid)
    if (.not.valid) return
    call euler_stage(candidate,temps,second_flux,2,second_diffusivity,valid)
    if (.not.valid) return
    do node=1,size(tree%nodes)
      if (owners(node)/=rank.or.tree%nodes(node)%level/=depth) cycle
      candidate%fields(node)%values=.5_dp*initial%fields(node)%values+.5_dp*candidate%fields(node)%values
      temps%fields(node)%values=first_temperature%fields(node)%values
      first_flux(node)%x=first_flux(node)%x+second_flux(node)%x
      first_flux(node)%y=first_flux(node)%y+second_flux(node)%y
      first_flux(node)%z=first_flux(node)%z+second_flux(node)%z
      if (.not.all(ieee_is_finite(first_flux(node)%x)).or..not.all(ieee_is_finite(first_flux(node)%y)).or. &
          .not.all(ieee_is_finite(first_flux(node)%z))) valid=.false.
    end do
    call agree(valid)
    if (.not.valid) return
    call recover_tree(candidate,temps,valid)
    if (.not.valid) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(temps%fields,temperature%fields)
    call move_alloc(first_flux,integrated_flux)
    maximum_diffusivity=max(first_diffusivity,second_diffusivity);ok=.true.
  contains
    subroutine recover_tree(u,t,success)
      type(amr_tree_state_3d), intent(in) :: u
      type(amr_tree_state_3d), intent(inout) :: t
      logical, intent(out) :: success
      call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,u,t,comm,success,depth=depth)
    end subroutine
    subroutine euler_stage(u,t,flux,stage,diffusivity,success)
      type(amr_tree_state_3d), intent(inout) :: u,t
      real(dp), intent(out) :: diffusivity
      integer, intent(in) :: stage
      real(dp) :: work_stamp
      logical, intent(out) :: success
      type(amr_tree_state_3d) :: packed,updated,ghosts
      type(amr_tree_flux_3d), allocatable, intent(out) :: flux(:)
      real(dp), allocatable :: fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:),fractions(:),les_rhs(:,:,:,:)
      real(dp) :: rhs(nvar),local_diffusivity,global_diffusivity
      character(len=32) :: local_faces(6)
      integer :: peer,n(3),dims(3),i,j,k,allocation_status,axis,level
      success=.false.;diffusivity=0
      if (present(fractions_start)) then
        if (stage==1) then
          fractions=fractions_start
        else
          fractions=fractions_end
        end if
      end if
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
        if (tree%nodes(peer)%level<depth) packed%fields(peer)%values=coarse_start%fields(peer)%values
      end do
      success=allocation_status==0;call agree(success)
      if (.not.success) return
      call copy_local(u,updated,success)
      if (.not.success) return
      call fill_mpi_amr_tree_ghosts_3d(tree,owners,packed,nvar+1,faces,ghosts,comm,success, &
        viscous_wall=.true.,inflow_state=inlet,wall_velocity=wall_velocity, &
        coarse_limiter=coarse_limiter,species_count=size(species), &
        coarse_future=coarse_end,coarse_fractions=fractions,depth=depth,ghost_width=width, &
        species=species,inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure)
      if (.not.success) return
      do peer=1,size(tree%nodes)
        if (owners(peer)/=rank.or.tree%nodes(peer)%level/=depth) cycle
        n=tree%nodes(peer)%shape;dims=n+4
        work_stamp=tree_work_start()
        associate(ghost=>ghosts%fields(peer)%values)
          call apply_amr_tree_pressure_ghosts_3d(tree,species,peer,faces,pressure,ghost,success,ghost_width=width)
          if (.not.success) exit
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
          allocate(les_rhs(nvar,n(1),n(2),n(3)),stat=allocation_status)
          if(allocation_status/=0) then
            success=.false.;exit
          end if
          if(dynamic) then
            call compute_dynamic_les_rhs_ghosted_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:), &
              h(:,peer),filter,filter_ratio,width,les_rhs,fx,fy,fz,success,grid,harmonic,local_diffusivity, &
              boundary_faces=local_faces,wall_temperature=wall_temperature)
          else
            call compute_les_ghosted_rhs_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:),n,h(:,peer), &
              cs,prt,sct,les_rhs,local_diffusivity,success,boundary_faces=local_faces, &
              wall_temperature=wall_temperature,model=selected_model,ci=isotropic,fx=fx,fy=fy,fz=fz)
          end if
          if (.not.success) exit
          success=dt*local_diffusivity*sum(1._dp/h(:,peer)**2)<=.2_dp
          if(.not.success) exit
          diffusivity=max(diffusivity,local_diffusivity)
          deallocate(les_rhs)
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
          allocate(flux(peer)%x,source=.5_dp*dt*fx,stat=allocation_status)
          if (allocation_status/=0) then
            success=.false.;exit
          end if
          allocate(flux(peer)%y,source=.5_dp*dt*fy,stat=allocation_status)
          if (allocation_status/=0) then
            success=.false.;exit
          end if
          allocate(flux(peer)%z,source=.5_dp*dt*fz,stat=allocation_status)
          if (allocation_status/=0) then
            success=.false.;exit
          end if
          deallocate(fx,fy,fz)
        end associate
        call tree_work_stop(peer,work_stamp)
      end do
      call agree(success)
      if (.not.success) return
      call MPI_Allreduce(diffusivity,global_diffusivity,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if (status/=MPI_SUCCESS) then
        success=.false.;return
      end if
      diffusivity=global_diffusivity
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
  end subroutine advance_mpi_amr_tree_les_level_3d
end module
