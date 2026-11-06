module les_mpi_amr_tree_mod
  use mpi_amr_tree_work_mod, only: tree_work_start,tree_work_stop
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use reactive_mol_3d_mod, only: convert_regular_energy_3d
  use les_amr_patch_tree_mod, only: convert_amr_les_tree_energy_3d
  use les_dynamic_smagorinsky_mod, only: compute_dynamic_les_rhs_ghosted_3d
  use les_filter_3d_mod, only: les_filter_weights
  use les_smagorinsky_3d_mod, only: compute_les_ghosted_rhs_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use amr_patch_tree_ghost_3d_mod, only: apply_amr_tree_pressure_ghosts_3d
  use mpi_amr_patch_tree_ghost_3d_mod, only: fill_mpi_amr_tree_ghosts_3d
  use mpi_amr_patch_tree_reflux_3d_mod, only: reflux_mpi_amr_tree_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: recover_mpi_amr_tree_temperatures_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  implicit none
  private
  public :: advance_mpi_amr_tree_les_3d,compute_mpi_amr_tree_les_timestep_3d,compute_mpi_amr_tree_les_rhs_3d
  public :: agree_les_frame_3d
contains
  subroutine compute_mpi_amr_tree_les_rhs_3d(tree,owners,species,state,temperature,root_spacing,faces, &
      cs,prt,sct,rates,flux,maximum_diffusivity,comm,ok, &
      inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter,outlet_pressure,model,ci, &
      test_model,test_ratio,grid_ratio,harmonic_mean,inflow_mach,inflow_pressure,target_depth,coarse_start,coarse_end,fractions, &
      frame,domain_lower,internal_energy)
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy
    character(len=:), allocatable :: sampling_limiter
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: root_spacing(3)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: cs,prt,sct
    real(dp), optional, intent(in) :: ci
    character(len=*), optional, intent(in) :: model
    real(dp) :: isotropic
    character(len=16) :: selected_model
    character(len=*), optional, intent(in) :: test_model
    integer, optional, intent(in) :: test_ratio,grid_ratio
    logical, optional, intent(in) :: harmonic_mean
    character(len=64) :: filter
    integer :: filter_ratio,grid,width
    logical :: dynamic,harmonic,filter_valid
    real(dp), allocatable :: weights(:)
    real(dp), intent(inout) :: maximum_diffusivity
    type(amr_tree_state_3d), intent(inout) :: rates
    type(amr_tree_flux_3d), allocatable, intent(inout) :: flux(:)
    type(amr_tree_state_3d), optional, intent(in) :: coarse_start,coarse_end
    real(dp), optional, intent(in) :: fractions(:)
    integer, optional, intent(in) :: target_depth
    integer :: depth
    integer :: presence(3),lowest(3),highest(3)
    type(amr_tree_flux_3d), allocatable :: result(:)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6)
    real(dp) :: pressure(6),mach(6),tin(6),pin(6)
    character(len=*), optional, intent(in) :: coarse_limiter
    type(amr_tree_state_3d) :: candidate,packed,ghosts,future_ordinary
    real(dp), allocatable :: h(:,:),inlet(:,:)
    real(dp) :: largest,diffusivity,global_diffusivity,lower(3),rotation(7),rotation_min(7),rotation_max(7)
    integer :: n(3),axis,level
    character(len=32) :: local_faces(6)
    real(dp), allocatable :: sampled_h(:,:,:)
    integer :: np,held_face,held_flag,held_low,held_high
    integer :: node,nvar,rank,status,ios,parent
    real(dp) :: controls(36),low(36),high(36)
    real(dp), allocatable :: reference_inlet(:,:)
    character(len=64) :: options(9),reference_options(9)
    integer :: flags(5),min_flags(5),max_flags(5),axis_min,axis_max
    real(dp) :: walls(24),reference_walls(24)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    logical :: valid
    if(present(coarse_limiter)) sampling_limiter=coarse_limiter
    if(present(internal_energy)) sampling_limiter='linear'
    ok=.false.
    nvar=reactive_nvar(size(species))
    call validate_mpi_amr_tree_local_3d(tree,owners,state,nvar,comm,valid)
    if (.not.valid) return
    call validate_mpi_amr_tree_local_3d(tree,owners,temperature,1,comm,valid)
    if (.not.valid) return
    held_flag=merge(1,0,present(internal_energy))
    call MPI_Allreduce(held_flag,held_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(held_flag,held_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.held_low/=held_high) return
    if(present(internal_energy)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,internal_energy,1,comm,valid)
      if(.not.valid) return
    endif
    np=nvar+1+held_flag
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
    depth=-1
    if(present(target_depth)) depth=target_depth
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
    valid=valid.and.(present(frame).eqv.present(domain_lower))
    if(present(frame)) then
      valid=valid.and.frame%is_valid()
      if(present(domain_lower)) valid=valid.and.all(ieee_is_finite(domain_lower))
    endif
    call agree(valid)
    if(.not.valid) return
    controls=[root_spacing,.2_dp,pressure,cs,prt,sct,isotropic,real(filter_ratio,dp),real(grid,dp), &
      merge(1._dp,0._dp,harmonic),real(depth,dp),mach,tin,pin]
    valid=all(ieee_is_finite(controls)).and.all(controls(:4)>0)
    valid=valid.and.all(pressure>=0).and.all(pressure<=1e15_dp)
    do node=1,6
      if (pressure(node)>0.and.faces(node)/='outflow') valid=.false.
    end do
    valid=valid.and.valid_reactive_boundaries_3d(faces)
    valid=valid.and.depth>=-1.and.depth<tree%level_count().and.cs>=0.and.cs<=2.and.min(prt,sct)>=1e-6_dp.and.max(prt,sct)<=1e6_dp
    valid=valid.and.isotropic>=0.and.isotropic<=1
    valid=valid.and.filter_valid.and. &
      (selected_model=='smagorinsky'.or.selected_model=='wale'.or.selected_model=='vreman'.or.dynamic)
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
    call MPI_Allreduce(controls,low,36,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,36,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(low/=high)) return
    flags=[merge(1,0,present(wall_temperature)),merge(1,0,present(wall_velocity)), &
      merge(1,0,present(coarse_limiter)),merge(1,0,present(frame)),merge(1,0,present(domain_lower))]
    call MPI_Allreduce(flags,min_flags,size(flags),MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,max_flags,size(flags),MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(min_flags/=max_flags)) return
    if(present(frame)) then
      call MPI_Allreduce(frame%axis,axis_min,1,MPI_INTEGER,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(frame%axis,axis_max,1,MPI_INTEGER,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.axis_min/=axis_max) return
      rotation=[domain_lower,frame%angular_speed,frame%axis_point]
      call MPI_Allreduce(rotation,rotation_min,size(rotation),MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(rotation,rotation_max,size(rotation),MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.any(rotation_min/=rotation_max)) return
    endif
    options(:6)=faces;options(7)='pcm';options(8)=selected_model;options(9)=filter
    if (present(coarse_limiter)) options(7)=coarse_limiter
    reference_options=options
    call MPI_Bcast(reference_options,64*9,MPI_CHARACTER,0,comm,status)
    if (status/=MPI_SUCCESS) return
    valid=all(options==reference_options)
    call agree(valid)
    if (.not.valid) return
    if (any(faces=='inflow'.and.mach<0)) then
      valid=present(inflow_state).and.present(inflow_temperature)
      if (valid) then
        valid=all(shape(inflow_state)==[nvar,6]).and.size(inflow_temperature)==6
        if (valid) valid=all(ieee_is_finite(inflow_state)).and.all(ieee_is_finite(inflow_temperature))
      end if
      call agree(valid)
      if (.not.valid) return
      allocate(inlet(np,6),reference_inlet(np,6),stat=ios)
      valid=ios==0;call agree(valid)
      if (.not.valid) return
      inlet(:nvar,:)=inflow_state;inlet(nvar+1,:)=inflow_temperature
      if(present(internal_energy)) then
        inlet(np,:)=0
        do held_face=1,6
          if(inlet(1,held_face)<=0) cycle
          inlet(np,held_face)=inlet(5,held_face)-.5_dp*sum(inlet(2:4,held_face)**2)/inlet(1,held_face)
        enddo
      endif
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
    presence=[merge(1,0,present(coarse_start)),merge(1,0,present(coarse_end)),merge(1,0,present(fractions))]
    call MPI_Allreduce(presence,lowest,3,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(presence,highest,3,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(lowest/=highest).or.any(presence/=presence(1))) return
    if(depth>0.and.presence(1)==0) return
    if(present(coarse_start)) then
      valid=depth>0.and.size(fractions)==depth
      if(valid) valid=all(ieee_is_finite(fractions)).and.all(fractions>=0).and.all(fractions<=1)
      call agree(valid)
      if(.not.valid) return
      call validate_mpi_amr_tree_local_3d(tree,owners,coarse_start,np,comm,valid,maximum_depth=depth-1)
      if(.not.valid) return
      call validate_mpi_amr_tree_local_3d(tree,owners,coarse_end,np,comm,valid,maximum_depth=depth-1)
      if(.not.valid) return
    end if
    allocate(h(3,size(tree%nodes)),packed%fields(size(tree%nodes)),result(size(tree%nodes)),stat=ios)
    valid=ios==0;call agree(valid)
    if(.not.valid) return
    h(:,1)=root_spacing
    do node=2,size(tree%nodes)
      parent=tree%nodes(node)%parent
      h(:,node)=h(:,parent)/real(tree%nodes(node)%patch%refinement_ratio,dp)
      if(any(h(:,node)<=0)) return
    end do
    candidate=state
    do node=1,size(tree%nodes)
      if(owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      allocate(packed%fields(node)%values(np,n(1),n(2),n(3)),stat=ios)
      if(ios/=0) exit
      packed%fields(node)%values(:nvar,:,:,:)=state%fields(node)%values
      packed%fields(node)%values(nvar+1,:,:,:)=temperature%fields(node)%values(1,:,:,:)
      if(present(internal_energy)) packed%fields(node)%values(np:np,:,:,:)=internal_energy%fields(node)%values
      if(tree%nodes(node)%level<depth) packed%fields(node)%values=coarse_start%fields(node)%values
      if(present(frame)) then
        ! Convert owner donors before the halo exchange and coarse interpolation.
        lower=domain_lower+real(origin(:,node),dp)*h(:,node)
        call convert_regular_energy_3d(frame,lower,h(:,node),packed%fields(node)%values,0,.false.,valid)
        if(.not.valid) exit
      endif
      candidate%fields(node)%values=0
    end do
    valid=valid.and.ios==0;call agree(valid)
    if(.not.valid) return
    if(present(frame).and.present(coarse_end)) then
      future_ordinary=coarse_end
      do node=1,size(tree%nodes)
        if(owners(node)/=rank.or.tree%nodes(node)%level>=depth) cycle
        lower=domain_lower+real(origin(:,node),dp)*h(:,node)
        call convert_regular_energy_3d(frame,lower,h(:,node),future_ordinary%fields(node)%values,0,.false.,valid)
        if(.not.valid) exit
      enddo
      call agree(valid)
      if(.not.valid) return
      call fill_mpi_amr_tree_ghosts_3d(tree,owners,packed,np,faces,ghosts,comm,valid, &
        viscous_wall=.true.,inflow_state=inlet,wall_velocity=wall_velocity, &
        coarse_limiter=sampling_limiter,species_count=size(species),ghost_width=width, &
        coarse_future=future_ordinary,coarse_fractions=fractions,depth=depth, &
        species=species,inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        stored_internal_energy=present(internal_energy),limit_internal_energy=present(internal_energy))
    else
      call fill_mpi_amr_tree_ghosts_3d(tree,owners,packed,np,faces,ghosts,comm,valid, &
        viscous_wall=.true.,inflow_state=inlet,wall_velocity=wall_velocity, &
        coarse_limiter=sampling_limiter,species_count=size(species),ghost_width=width, &
        coarse_future=coarse_end,coarse_fractions=fractions,depth=depth, &
        species=species,inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        stored_internal_energy=present(internal_energy),limit_internal_energy=present(internal_energy))
    endif
    if(.not.valid) return
    largest=0
    do node=1,size(tree%nodes)
      if(owners(node)/=rank) cycle
      if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
      n=tree%nodes(node)%shape
      associate(ghost=>ghosts%fields(node)%values)
        call apply_amr_tree_pressure_ghosts_3d(tree,species,node,faces,pressure,ghost,valid,ghost_width=width, &
        stored_internal_energy=present(internal_energy))
        if(.not.valid) exit
        local_faces='outflow';level=tree%nodes(node)%level+1
        do axis=1,3
          if(origin(axis,node)==0) local_faces(2*axis-1)=faces(2*axis-1)
          if(origin(axis,node)+int(n(axis),int64)==domain(axis,level)) local_faces(2*axis)=faces(2*axis)
        end do
        allocate(result(node)%x(nvar,n(1)+1,n(2),n(3)),result(node)%y(nvar,n(1),n(2)+1,n(3)), &
          result(node)%z(nvar,n(1),n(2),n(3)+1),stat=ios)
        if(ios/=0) then
          valid=.false.;exit
        end if
      if(present(internal_energy)) sampled_h=ghost(np,:,:,:)
        if(dynamic) then
          call compute_dynamic_les_rhs_ghosted_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:), &
            h(:,node),filter,filter_ratio,width,candidate%fields(node)%values, &
            result(node)%x,result(node)%y,result(node)%z,valid,grid,harmonic,diffusivity, &
            boundary_faces=local_faces,wall_temperature=wall_temperature, &
              internal_energy_density=sampled_h)
        else
          call compute_les_ghosted_rhs_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:),n,h(:,node), &
            cs,prt,sct,candidate%fields(node)%values,diffusivity,valid,boundary_faces=local_faces, &
            wall_temperature=wall_temperature,model=selected_model,ci=isotropic, &
            fx=result(node)%x,fy=result(node)%y,fz=result(node)%z,internal_energy_density=sampled_h)
        end if
        if(.not.valid) exit
        largest=max(largest,diffusivity)
      end associate
    end do
    call agree(valid)
    if(.not.valid) return
    call MPI_Allreduce(largest,global_diffusivity,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    call move_alloc(candidate%fields,rates%fields)
    call move_alloc(result,flux)
    maximum_diffusivity=global_diffusivity;ok=.true.
  contains
    subroutine agree(condition)
      logical, intent(inout) :: condition
      logical :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      condition=.false.
      if(ierr==MPI_SUCCESS) condition=result
    end subroutine agree
  end subroutine compute_mpi_amr_tree_les_rhs_3d
  recursive subroutine compute_mpi_amr_tree_les_timestep_3d(tree,owners,species,state,temperature,root_spacing,cfl,faces, &
      cs,prt,sct,dt,maximum_diffusivity,comm,ok, &
      inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter,outlet_pressure,model,ci, &
      test_model,test_ratio,grid_ratio,harmonic_mean,subcycle,inflow_mach,inflow_pressure,frame,domain_lower,internal_energy)
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy
    character(len=:), allocatable :: sampling_limiter
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: root_spacing(3),cfl
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: cs,prt,sct
    real(dp), optional, intent(in) :: ci
    character(len=*), optional, intent(in) :: model
    real(dp) :: isotropic
    character(len=16) :: selected_model
    character(len=*), optional, intent(in) :: test_model
    integer, optional, intent(in) :: test_ratio,grid_ratio
    logical, optional, intent(in) :: harmonic_mean,subcycle
    character(len=64) :: filter
    integer :: filter_ratio,grid,width
    logical :: dynamic,harmonic,filter_valid,subcycling
    real(dp), allocatable :: weights(:)
    real(dp), intent(inout) :: dt,maximum_diffusivity
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6)
    real(dp) :: pressure(6),mach(6),tin(6),pin(6)
    character(len=*), optional, intent(in) :: coarse_limiter
    type(amr_tree_state_3d), allocatable :: timestep_h
    type(amr_tree_state_3d) :: initial,temps,packed,ghosts,ordinary
    real(dp), allocatable :: h(:,:),inlet(:,:)
    real(dp) :: candidate_dt,largest,diffusivity,rate,global_dt,global_diffusivity
    real(dp), allocatable :: rhs(:,:,:,:)
    integer :: n(3),axis,level
    character(len=32) :: local_faces(6)
    real(dp), allocatable :: sampled_h(:,:,:)
    integer :: np,held_face,held_flag,held_low,held_high
    integer :: node,nvar,rank,status,ios,parent
    real(dp) :: controls(36),low(36),high(36)
    real(dp), allocatable :: reference_inlet(:,:)
    character(len=64) :: options(9),reference_options(9)
    integer :: flags(3),min_flags(3),max_flags(3)
    real(dp) :: walls(24),reference_walls(24)
    integer(int64), allocatable :: origin(:,:),domain(:,:),counts(:)
    integer, allocatable :: ratios(:)
    integer(int64) :: factor
    logical :: valid
    if(present(coarse_limiter)) sampling_limiter=coarse_limiter
    if(present(internal_energy)) sampling_limiter='linear'
    ok=.false.
    nvar=reactive_nvar(size(species))
    held_flag=merge(1,0,present(internal_energy))
    call MPI_Allreduce(held_flag,held_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(held_flag,held_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.held_low/=held_high) return
    if(present(internal_energy)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,internal_energy,1,comm,valid)
      if(.not.valid) return
    endif
    np=nvar+1+held_flag
    call agree_les_frame_3d(frame,domain_lower,comm,valid)
    if(.not.valid) return
    if(present(frame)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,state,nvar,comm,valid)
      if(.not.valid) return
      call MPI_Comm_rank(comm,rank,status)
      if(status/=MPI_SUCCESS) return
      ordinary=state
      call convert_amr_les_tree_energy_3d(tree,ordinary,root_spacing,frame,domain_lower,.false.,valid,owners,rank)
      call agree(valid)
      if(.not.valid) return
      call compute_mpi_amr_tree_les_timestep_3d(tree,owners,species,ordinary,temperature,root_spacing,cfl,faces, &
        cs,prt,sct,dt,maximum_diffusivity,comm,ok, &
        inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter,outlet_pressure,model,ci, &
        test_model,test_ratio,grid_ratio,harmonic_mean,subcycle,inflow_mach,inflow_pressure,internal_energy=internal_energy)
      return
    endif
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
    subcycling=.false.
    if(present(subcycle)) subcycling=subcycle
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
    controls=[root_spacing,cfl,pressure,cs,prt,sct,isotropic,real(filter_ratio,dp),real(grid,dp), &
      merge(1._dp,0._dp,harmonic),merge(1._dp,0._dp,subcycling),mach,tin,pin]
    valid=all(ieee_is_finite(controls)).and.all(controls(:4)>0)
    valid=valid.and.all(pressure>=0).and.all(pressure<=1e15_dp)
    do node=1,6
      if (pressure(node)>0.and.faces(node)/='outflow') valid=.false.
    end do
    valid=valid.and.valid_reactive_boundaries_3d(faces)
    valid=valid.and.cfl<=.2_dp.and.cs>=0.and.cs<=2.and.min(prt,sct)>=1e-6_dp.and.max(prt,sct)<=1e6_dp
    valid=valid.and.isotropic>=0.and.isotropic<=1
    valid=valid.and.filter_valid.and. &
      (selected_model=='smagorinsky'.or.selected_model=='wale'.or.selected_model=='vreman'.or.dynamic)
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
    call MPI_Allreduce(controls,low,36,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,36,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
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
    if (any(faces=='inflow'.and.mach<0)) then
      valid=present(inflow_state).and.present(inflow_temperature)
      if (valid) then
        valid=all(shape(inflow_state)==[nvar,6]).and.size(inflow_temperature)==6
        if (valid) valid=all(ieee_is_finite(inflow_state)).and.all(ieee_is_finite(inflow_temperature))
      end if
      call agree(valid)
      if (.not.valid) return
      allocate(inlet(np,6),reference_inlet(np,6),stat=ios)
      valid=ios==0;call agree(valid)
      if (.not.valid) return
      inlet(:nvar,:)=inflow_state;inlet(nvar+1,:)=inflow_temperature
      if(present(internal_energy)) then
        inlet(np,:)=0
        do held_face=1,6
          if(inlet(1,held_face)<=0) cycle
          inlet(np,held_face)=inlet(5,held_face)-.5_dp*sum(inlet(2:4,held_face)**2)/inlet(1,held_face)
        enddo
      endif
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
    allocate(counts(size(ratios)),stat=ios)
    valid=ios==0;call agree(valid)
    if(.not.valid) return
    counts=1
    if(subcycling) then
      do level=2,size(ratios)
        factor=int(ratios(level),int64)**2
        if(counts(level-1)>huge(1_int64)/factor) then
          valid=.false.;exit
        end if
        counts(level)=counts(level-1)*factor
      end do
    end if
    call agree(valid)
    if(.not.valid) return
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
    if(present(internal_energy)) then
      allocate(timestep_h);timestep_h=internal_energy
      call average_down_mpi_amr_tree_3d(tree,owners,timestep_h,1,comm,valid)
      if(.not.valid) return
    endif
    call recover_tree(initial,temps,valid)
    if (.not.valid) return
    allocate(packed%fields(size(tree%nodes)),stat=ios)
    valid=ios==0;call agree(valid)
    if (.not.valid) return
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      allocate(packed%fields(node)%values(np,n(1),n(2),n(3)),stat=ios)
      if (ios/=0) exit
      packed%fields(node)%values(:nvar,:,:,:)=initial%fields(node)%values
      packed%fields(node)%values(nvar+1,:,:,:)=temps%fields(node)%values(1,:,:,:)
      if(present(internal_energy)) packed%fields(node)%values(np:np,:,:,:)=timestep_h%fields(node)%values
    end do
    valid=ios==0;call agree(valid)
    if (.not.valid) return
    call fill_mpi_amr_tree_ghosts_3d(tree,owners,packed,np,faces,ghosts,comm,valid, &
      viscous_wall=.true.,inflow_state=inlet,wall_velocity=wall_velocity, &
      coarse_limiter=sampling_limiter,species_count=size(species),ghost_width=width, &
        species=species,inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        stored_internal_energy=present(internal_energy),limit_internal_energy=present(internal_energy))
    if (.not.valid) return
    candidate_dt=1e6_dp;largest=0
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      associate(ghost=>ghosts%fields(node)%values)
        call apply_amr_tree_pressure_ghosts_3d(tree,species,node,faces,pressure,ghost,valid,ghost_width=width, &
        stored_internal_energy=present(internal_energy))
        if (.not.valid) exit
        local_faces='outflow';level=tree%nodes(node)%level+1
        do axis=1,3
          if (origin(axis,node)==0) local_faces(2*axis-1)=faces(2*axis-1)
          if (origin(axis,node)+int(n(axis),int64)==domain(axis,level)) local_faces(2*axis)=faces(2*axis)
        end do
        allocate(rhs(nvar,n(1),n(2),n(3)),stat=ios)
        if (ios/=0) then
          valid=.false.;exit
        end if
      if(present(internal_energy)) sampled_h=ghost(np,:,:,:)
        if(dynamic) then
          block
            real(dp), allocatable :: fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:)
            allocate(fx(nvar,n(1)+1,n(2),n(3)),fy(nvar,n(1),n(2)+1,n(3)), &
              fz(nvar,n(1),n(2),n(3)+1),stat=ios)
            valid=ios==0
            if(valid) call compute_dynamic_les_rhs_ghosted_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:), &
              h(:,node),filter,filter_ratio,width,rhs,fx,fy,fz,valid,grid,harmonic,diffusivity, &
              boundary_faces=local_faces,wall_temperature=wall_temperature, &
              internal_energy_density=sampled_h)
          end block
        else
        call compute_les_ghosted_rhs_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:),n,h(:,node), &
          cs,prt,sct,rhs,diffusivity,valid,boundary_faces=local_faces,wall_temperature=wall_temperature, &
          model=selected_model,ci=isotropic,internal_energy_density=sampled_h)
        end if
        if (.not.valid) exit
        rate=diffusivity*sum(1._dp/h(:,node)**2)/real(counts(level),dp)
        if (rate>cfl/candidate_dt) candidate_dt=cfl/rate
        largest=max(largest,diffusivity)
        deallocate(rhs)
      end associate
    end do
    call agree(valid)
    if (.not.valid) return
    call MPI_Allreduce(candidate_dt,global_dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(largest,global_diffusivity,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    dt=global_dt;maximum_diffusivity=global_diffusivity;ok=.true.
  contains
    subroutine recover_tree(u,t,success)
      type(amr_tree_state_3d), intent(in) :: u
      type(amr_tree_state_3d), intent(inout) :: t
      logical, intent(out) :: success
      call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,u,t,comm,success,internal_energy=timestep_h)
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
  end subroutine compute_mpi_amr_tree_les_timestep_3d
  recursive subroutine advance_mpi_amr_tree_les_3d(tree,owners,species,state,temperature,root_spacing,dt,faces, &
      cs,prt,sct,boundary_change,maximum_correction,comm,ok, &
      inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter,outlet_pressure,model,ci, &
      test_model,test_ratio,grid_ratio,harmonic_mean,inflow_mach,inflow_pressure,frame,domain_lower)
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: root_spacing(3),dt
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: cs,prt,sct
    real(dp), optional, intent(in) :: ci
    character(len=*), optional, intent(in) :: model
    real(dp) :: isotropic
    character(len=16) :: selected_model
    character(len=*), optional, intent(in) :: test_model
    integer, optional, intent(in) :: test_ratio,grid_ratio
    logical, optional, intent(in) :: harmonic_mean
    character(len=64) :: filter
    integer :: filter_ratio,grid,width
    logical :: dynamic,harmonic,filter_valid
    real(dp), allocatable :: weights(:)
    real(dp), intent(inout) :: boundary_change(:),maximum_correction
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6)
    real(dp) :: pressure(6),mach(6),tin(6),pin(6)
    character(len=*), optional, intent(in) :: coarse_limiter
    type(amr_tree_state_3d) :: initial,candidate,temps,first_temperature,ordinary,ordinary_t
    real(dp), allocatable :: h(:,:),inlet(:,:)
    real(dp) :: first_budget(size(boundary_change)),second_budget(size(boundary_change)), &
      first_max,second_max
    integer :: node,nvar,rank,status,ios,parent
    real(dp) :: controls(35),low(35),high(35)
    real(dp), allocatable :: reference_inlet(:,:)
    character(len=64) :: options(9),reference_options(9)
    integer :: flags(3),min_flags(3),max_flags(3)
    real(dp) :: walls(24),reference_walls(24)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    logical :: valid
    real(dp) :: rotated_budget(size(boundary_change)),rotated_max
    ok=.false.
    nvar=reactive_nvar(size(species))
    call agree_les_frame_3d(frame,domain_lower,comm,valid)
    if(.not.valid) return
    if(present(frame)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,state,nvar,comm,valid)
      if(.not.valid) return
      call MPI_Comm_rank(comm,rank,status)
      if(status/=MPI_SUCCESS) return
      ordinary=state;ordinary_t=temperature
      call convert_amr_les_tree_energy_3d(tree,ordinary,root_spacing,frame,domain_lower,.false.,valid,owners,rank)
      call agree(valid)
      if(.not.valid) return
      call advance_mpi_amr_tree_les_3d(tree,owners,species,ordinary,ordinary_t,root_spacing,dt,faces, &
        cs,prt,sct,rotated_budget,rotated_max,comm,valid, &
        inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter,outlet_pressure,model,ci, &
        test_model,test_ratio,grid_ratio,harmonic_mean,inflow_mach,inflow_pressure)
      if(.not.valid) return
      call convert_amr_les_tree_energy_3d(tree,ordinary,root_spacing,frame,domain_lower,.true.,valid,owners,rank)
      call agree(valid)
      if(.not.valid) return
      state=ordinary;temperature=ordinary_t;boundary_change=rotated_budget;maximum_correction=rotated_max;ok=.true.
      return
    endif
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
    controls=[root_spacing,dt,pressure,cs,prt,sct,isotropic,real(filter_ratio,dp),real(grid,dp),merge(1._dp, &
      0._dp,harmonic),mach,tin,pin]
    valid=all(ieee_is_finite(controls)).and.all(controls(:4)>0)
    valid=valid.and.all(pressure>=0).and.all(pressure<=1e15_dp)
    do node=1,6
      if (pressure(node)>0.and.faces(node)/='outflow') valid=.false.
    end do
    valid=valid.and.size(boundary_change)==nvar.and.valid_reactive_boundaries_3d(faces)
    valid=valid.and.dt<=1e6_dp.and.cs>=0.and.cs<=2.and.min(prt,sct)>=1e-6_dp.and.max(prt,sct)<=1e6_dp
    valid=valid.and.isotropic>=0.and.isotropic<=1
    valid=valid.and.filter_valid.and. &
      (selected_model=='smagorinsky'.or.selected_model=='wale'.or.selected_model=='vreman'.or.dynamic)
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
    call average_down_mpi_amr_tree_3d(tree,owners,initial,nvar,comm,valid)
    if (.not.valid) return
    call recover_tree(initial,temps,valid)
    if (.not.valid) return
    call copy_local(initial,candidate,valid)
    if (.not.valid) return
    call euler_stage(candidate,temps,first_budget,first_max,valid)
    if (.not.valid) return
    call copy_local(temps,first_temperature,valid)
    if (.not.valid) return
    call euler_stage(candidate,temps,second_budget,second_max,valid)
    if (.not.valid) return
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      candidate%fields(node)%values=.5_dp*initial%fields(node)%values+.5_dp*candidate%fields(node)%values
    end do
    call average_down_mpi_amr_tree_3d(tree,owners,candidate,nvar,comm,valid)
    if (.not.valid) return
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      temps%fields(node)%values=first_temperature%fields(node)%values
    end do
    call recover_tree(candidate,temps,valid)
    if (.not.valid) return
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(temps%fields,temperature%fields)
    boundary_change=.5_dp*(first_budget+second_budget)
    maximum_correction=max(first_max,second_max);ok=.true.
  contains
    subroutine recover_tree(u,t,success)
      type(amr_tree_state_3d), intent(in) :: u
      type(amr_tree_state_3d), intent(inout) :: t
      logical, intent(out) :: success
      call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,u,t,comm,success)
    end subroutine
    subroutine euler_stage(u,t,budget,largest,success)
      type(amr_tree_state_3d), intent(inout) :: u,t
      real(dp), intent(out) :: budget(:),largest
      real(dp) :: work_stamp
      logical, intent(out) :: success
      type(amr_tree_state_3d) :: packed,updated,ghosts
      type(amr_tree_flux_3d), allocatable :: flux(:)
      real(dp), allocatable :: fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:),les_rhs(:,:,:,:)
      real(dp) :: correction(nvar),rhs(nvar),diffusivity
      character(len=32) :: local_faces(6)
      integer :: peer,n(3),dims(3),i,j,k,allocation_status,axis,level
      success=.false.;budget=0;largest=0;correction=0
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
        coarse_limiter=coarse_limiter,species_count=size(species),ghost_width=width, &
        species=species,inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure)
      if (.not.success) return
      do peer=1,size(tree%nodes)
        if (owners(peer)/=rank) cycle
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
          if (allocation_status/=0) then
            success=.false.;exit
          end if
          if(dynamic) then
            call compute_dynamic_les_rhs_ghosted_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:), &
              h(:,peer),filter,filter_ratio,width,les_rhs,fx,fy,fz,success,grid,harmonic,diffusivity, &
              boundary_faces=local_faces,wall_temperature=wall_temperature)
          else
          call compute_les_ghosted_rhs_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:),n,h(:,peer), &
            cs,prt,sct,les_rhs,diffusivity,success,boundary_faces=local_faces, &
            wall_temperature=wall_temperature,model=selected_model,ci=isotropic,fx=fx,fy=fy,fz=fz)
          end if
          if (.not.success) exit
          success=dt*diffusivity*sum(1._dp/h(:,peer)**2)<=.2_dp
          if (.not.success) exit
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
        call tree_work_stop(peer,work_stamp)
      end do
      call agree(success)
      if (.not.success) return
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
  end subroutine advance_mpi_amr_tree_les_3d
  subroutine agree_les_frame_3d(frame,domain_lower,comm,ok)
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer :: flags(2),low_flags(2),high_flags(2),axis_min,axis_max,status
    real(dp) :: values(7),minimum(7),maximum(7)
    logical :: valid,agreed
    ok=.false.
    flags=[merge(1,0,present(frame)),merge(1,0,present(domain_lower))]
    call MPI_Allreduce(flags,low_flags,2,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high_flags,2,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(low_flags/=high_flags).or.flags(1)/=flags(2)) return
    if(.not.present(frame)) then
      ok=.true.;return
    endif
    valid=frame%is_valid().and.all(ieee_is_finite(domain_lower))
    call MPI_Allreduce(valid,agreed,1,MPI_LOGICAL,MPI_LAND,comm,status)
    if(status/=MPI_SUCCESS) return
    if(.not.agreed) return
    call MPI_Allreduce(frame%axis,axis_min,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(frame%axis,axis_max,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(axis_min/=axis_max) return
    values=[frame%angular_speed,frame%axis_point,domain_lower]
    call MPI_Allreduce(values,minimum,7,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(values,maximum,7,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(minimum/=maximum)) return
    ok=.true.
  end subroutine agree_les_frame_3d
end module
