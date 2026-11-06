module mpi_amr_patch_tree_transport_3d_mod
  use mpi_amr_tree_work_mod, only: tree_work_start,tree_work_stop
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use gas_transport_mod, only: gas_transport_species
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,nasa7_eos_bounds
  use reactive_1d_mod, only: reactive_nvar
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use reactive_mol_3d_mod, only: convert_regular_energy_3d
  use reactive_transport_3d_mod, only: physical_transport_ghosted_fluxes_3d, &
    reactive_transport_timestep_3d,reactive_isothermal_timestep_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use amr_patch_tree_ghost_3d_mod, only: apply_amr_tree_pressure_ghosts_3d
  use mpi_amr_patch_tree_ghost_3d_mod, only: fill_mpi_amr_tree_ghosts_3d
  use mpi_amr_patch_tree_reflux_3d_mod, only: reflux_mpi_amr_tree_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: recover_mpi_amr_tree_temperatures_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  use reactive_transport_3d_mod, only: valid_constant_transport_coefficients
  implicit none
  private
  public :: advance_mpi_amr_tree_transport_3d,compute_mpi_amr_tree_transport_timestep_3d
  public :: compute_sdc_mpi_amr_tree_transport_3d
  public :: agree_mpi_constant_transport_coefficients
  public :: agree_mpi_thermal_slip_faces
contains
  subroutine agree_mpi_thermal_slip_faces(requested,comm,effective,ok,faces,conduction,wall_temperature)
    logical, optional, intent(in) :: requested(6)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: effective(6),ok
    character(len=*), optional, intent(in) :: faces(6)
    logical, intent(in) :: conduction
    real(dp), optional, intent(in) :: wall_temperature(:)
    integer :: values(7),low(7),high(7),status
    logical :: valid
    ok=.false.;effective=.false.;values=0;valid=.true.
    if(present(requested)) values(:6)=merge(1,0,requested)
    if(any(values(:6)==1)) then
      if(.not.present(faces).or..not.conduction.or..not.present(wall_temperature)) valid=.false.
      if(present(faces)) then
        if(any(values(:6)==1.and.faces/='slip_wall')) valid=.false.
      endif
      if(present(wall_temperature)) then
        if(size(wall_temperature)/=6) then
          valid=.false.
        else
          if(any(values(:6)==1.and.(.not.ieee_is_finite(wall_temperature).or.wall_temperature<=0))) valid=.false.
        endif
      endif
    endif
    values(7)=merge(1,0,valid)
    call MPI_Allreduce(values,low,7,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(values,high,7,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low(:6)/=high(:6)).or.low(7)==0) return
    effective=values(:6)==1;ok=.true.
  end subroutine agree_mpi_thermal_slip_faces

  subroutine agree_mpi_constant_transport_coefficients(coefficients,comm,ok,species)
    real(dp), optional, intent(in) :: coefficients(:)
    type(nasa7_species), optional, intent(in) :: species(:)
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    integer :: flags(4),lo(4),hi(4),status
    real(dp) :: values(7),low(7),high(7)
    flags=0;values=0;ok=.false.
    if(present(coefficients)) then
      flags(1)=size(coefficients)
      if(valid_constant_transport_coefficients(coefficients)) then
        flags(2)=1;values(1:size(coefficients))=coefficients
      endif
    else
      flags(2)=1
    endif
    flags(4)=1
    if(present(species)) then
      flags(3)=1;flags(4)=0
      if(size(species)>0) then
        if(all(species%pinned_eos_extrapolation.eqv.species(1)%pinned_eos_extrapolation)) then
          flags(4)=1
          if(species(1)%pinned_eos_extrapolation) then
            values(5)=1
            block
              logical :: policy_ok
              call nasa7_eos_bounds(species,values(6),values(7),policy_ok)
              flags(4)=merge(1,0,policy_ok)
            end block
          endif
        endif
      endif
    endif
    call MPI_Allreduce(flags,lo,4,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,hi,4,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(lo/=hi).or.(lo(2)/=1.or.lo(4)/=1)) return
    call MPI_Allreduce(values,low,7,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(values,high,7,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.any(low/=high)) return
    ok=.true.
  end subroutine agree_mpi_constant_transport_coefficients
  subroutine compute_mpi_amr_tree_transport_timestep_3d(tree,owners,species,transport,state,temperature, &
      root_spacing,cfl,viscosity,conduction,diffusion,dt,maximum_diffusivity,comm,ok,faces,wall_temperature, &
      constant_coefficients,internal_energy,isothermal_slip_faces)
    logical, optional, intent(in) :: isothermal_slip_faces(6)
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy
    real(dp), optional, intent(in) :: constant_coefficients(:)
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
    real(dp), allocatable :: sampled_h(:,:,:)
    real(dp), allocatable :: h(:,:)
    character(len=*), optional, intent(in) :: faces(6)
    real(dp), optional, intent(in) :: wall_temperature(:)
    real(dp) :: candidate,local_dt,local_diffusivity,diffusivity,wall_dt
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    character(len=32) :: local_faces(6)
    integer :: held_flag,held_low,held_high,np,held_face
    integer :: node,n(3),axis,level,rank,status,ios,parent
    integer :: flags(5),min_flags(5),max_flags(5)
    real(dp) :: controls(10),low(10),high(10),global_dt,global_diffusivity
    character(len=64) :: options(6),reference_options(6)
    logical :: thermal_slip(6),local_thermal_slip(6)
    logical :: valid
    ok=.false.
    call agree_mpi_thermal_slip_faces(isothermal_slip_faces,comm,thermal_slip,ok,faces,conduction,wall_temperature)
    if(.not.ok) return
    ok=.false.
    call agree_mpi_constant_transport_coefficients(constant_coefficients,comm,ok,species=species)
    if(.not.ok) return
    ok=.false.
    call validate_mpi_amr_tree_local_3d(tree,owners,state,reactive_nvar(size(species)),comm,valid)
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
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    controls=0;controls(:4)=[root_spacing,cfl]
    valid=all(ieee_is_finite(controls)).and.all(controls(:4)>0).and.size(transport)==size(species)
    if(any(thermal_slip)) then
      if(.not.present(faces).or..not.conduction.or..not.present(wall_temperature)) valid=.false.
      if(present(faces)) then
        if(any(thermal_slip.and.faces/='slip_wall')) valid=.false.
      endif
    endif
    if (present(faces)) then
      if (.not.valid_reactive_boundaries_3d(faces)) valid=.false.
      if (any(faces=='isothermal_wall').or.any(thermal_slip)) then
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
      if(present(internal_energy)) sampled_h=internal_energy%fields(node)%values(1,:,:,:)
      call reactive_transport_timestep_3d(species,transport,state%fields(node)%values,temperature%fields(node)%values(1,:,:,:), &
        n(1),n(2),n(3),h(1,node),h(2,node),h(3,node),cfl,viscosity,conduction,diffusion,local_dt,local_diffusivity,valid, &
          constant_coefficients=constant_coefficients,internal_energy_density=sampled_h)
      if (.not.valid) exit
      candidate=min(candidate,local_dt);diffusivity=max(diffusivity,local_diffusivity)
      if (present(faces)) then
        if (any(faces=='isothermal_wall').or.any(thermal_slip)) then
          local_faces='outflow';level=tree%nodes(node)%level+1
          do axis=1,3
            if (origin(axis,node)==0.and. &
                (faces(2*axis-1)=='isothermal_wall'.or.thermal_slip(2*axis-1))) &
              local_faces(2*axis-1)=faces(2*axis-1)
            if (origin(axis,node)+int(n(axis),int64)==domain(axis,level).and. &
                (faces(2*axis)=='isothermal_wall'.or.thermal_slip(2*axis))) &
              local_faces(2*axis)=faces(2*axis)
          end do
          local_thermal_slip=thermal_slip.and.local_faces=='slip_wall'
          call reactive_isothermal_timestep_3d(species,transport,state%fields(node)%values, &
            temperature%fields(node)%values(1,:,:,:),h(:,node),local_faces,wall_temperature,cfl,wall_dt,valid, &
          constant_coefficients=constant_coefficients,internal_energy_density=sampled_h,isothermal_slip_faces=local_thermal_slip)
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
  subroutine compute_sdc_mpi_amr_tree_transport_3d(tree,owners,species,transport,state,temperature,root_spacing,dt,faces, &
      flags,transport_rate,integrated_flux,minimum_theta,comm,ok,inflow_state,inflow_temperature, &
      wall_temperature,wall_velocity,outlet_pressure,inflow_mach,inflow_pressure, &
      target_depth,coarse_start,coarse_end,coarse_fractions, &
      constant_coefficients,frame,domain_lower,internal_energy,total_energy_primitives,isothermal_slip_faces)
    logical, optional, intent(in) :: isothermal_slip_faces(6)
    logical, optional, intent(in) :: total_energy_primitives
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy
    real(dp), optional, intent(in) :: constant_coefficients(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: root_spacing(3),dt
    character(len=*), intent(in) :: faces(6)
    logical, intent(in) :: flags(4)
    type(amr_tree_state_3d), intent(inout) :: transport_rate
    type(amr_tree_flux_3d), allocatable, intent(inout) :: integrated_flux(:)
    real(dp), intent(inout) :: minimum_theta
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6)
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d), optional, intent(in) :: coarse_start,coarse_end
    real(dp), optional, intent(in) :: coarse_fractions(:)
    type(amr_tree_state_3d) :: local_state,local_temperature
    real(dp) :: budget(5+size(species)),maximum
    local_state=state;local_temperature=temperature
    call advance_mpi_amr_tree_transport_3d(tree,owners,species,transport,local_state,local_temperature,root_spacing,dt,faces, &
      flags(1),flags(2),flags(3),flags(4),budget,maximum,minimum_theta,comm,ok, &
      inflow_state,inflow_temperature,wall_temperature,wall_velocity,'mc',outlet_pressure,inflow_mach,inflow_pressure, &
      transport_rate,integrated_flux,target_depth,coarse_start,coarse_end,coarse_fractions, &
          constant_coefficients=constant_coefficients,frame=frame,domain_lower=domain_lower,internal_energy=internal_energy, &
          total_energy_primitives=total_energy_primitives,isothermal_slip_faces=isothermal_slip_faces)
  end subroutine compute_sdc_mpi_amr_tree_transport_3d

  subroutine advance_mpi_amr_tree_transport_3d(tree,owners,species,transport,state,temperature,root_spacing,dt,faces, &
      viscosity,conduction,diffusion,barodiffusion,boundary_change,maximum_correction,minimum_theta,comm,ok, &
      inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter,outlet_pressure,inflow_mach,inflow_pressure, &
      transport_rate,integrated_flux,target_depth,coarse_start,coarse_end,coarse_fractions, &
      constant_coefficients,frame,domain_lower,internal_energy,total_energy_primitives,isothermal_slip_faces)
    logical, optional, intent(in) :: isothermal_slip_faces(6)
    logical, optional, intent(in) :: total_energy_primitives
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy
    real(dp), optional, intent(in) :: constant_coefficients(:)
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3)
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
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6)
    type(amr_tree_state_3d), optional, intent(inout) :: transport_rate
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: integrated_flux(:)
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d), optional, intent(in) :: coarse_start,coarse_end
    real(dp), optional, intent(in) :: coarse_fractions(:)
    real(dp) :: pressure(6),mach(6),pin(6),tin(6)
    character(len=:), allocatable :: sampling_limiter
    character(len=*), optional, intent(in) :: coarse_limiter
    type(amr_tree_state_3d) :: initial,candidate,temps,initial_temperature
    real(dp), allocatable :: h(:,:),inlet(:,:)
    real(dp) :: first_budget(size(boundary_change)),second_budget(size(boundary_change)), &
      first_max,second_max,first_theta,second_theta
    integer :: held_flag,held_low,held_high,np,held_face
    integer :: node,nvar,rank,status,ios,parent,depth,selection(4),selection_low(4),selection_high(4)
    real(dp) :: controls(28),low(28),high(28)
    real(dp) :: rotation(7),rotation_low(7),rotation_high(7)
    integer :: axis_low,axis_high
    real(dp), allocatable :: reference_inlet(:,:)
    character(len=64) :: options(7),reference_options(7)
    integer :: flags(9),min_flags(9),max_flags(9)
    real(dp) :: walls(24),reference_walls(24)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    logical :: thermal_slip(6),local_thermal_slip(6)
    logical :: valid,stage_only,use_total_energy
    if(present(coarse_limiter)) sampling_limiter=coarse_limiter
    if(present(internal_energy)) sampling_limiter='linear'
    ok=.false.
    call agree_mpi_thermal_slip_faces(isothermal_slip_faces,comm,thermal_slip,ok,faces,conduction,wall_temperature)
    if(.not.ok) return
    ok=.false.
    call agree_mpi_constant_transport_coefficients(constant_coefficients,comm,ok,species=species)
    if(.not.ok) return
    ok=.false.
    use_total_energy=.false.
    if(present(total_energy_primitives)) use_total_energy=total_energy_primitives
    held_flag=merge(1,0,use_total_energy)
    call MPI_Allreduce(held_flag,held_low,1,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(held_flag,held_high,1,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS.or.held_low/=held_high) return
    stage_only=present(transport_rate)
    valid=present(frame).eqv.present(domain_lower)
    if(present(frame)) then
      valid=valid.and.stage_only.and.frame%is_valid()
      if(present(domain_lower)) valid=valid.and.all(ieee_is_finite(domain_lower))
    endif
    call agree(valid)
    if(.not.valid) return
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
    call MPI_Comm_rank(comm,rank,status)
    if (status/=MPI_SUCCESS) return
    np=nvar+1+held_flag
    depth=-1
    if(present(target_depth)) depth=target_depth
    selection=[depth,merge(1,0,present(coarse_start)),merge(1,0,present(coarse_end)), &
      merge(1,0,present(coarse_fractions))]
    call MPI_Allreduce(selection,selection_low,4,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(selection,selection_high,4,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(selection_low/=selection_high)) return
    valid=depth>=-1.and.depth<tree%level_count()
    if(present(target_depth)) valid=valid.and.depth>=0.and.stage_only
    valid=valid.and.all(selection(2:4)==selection(2))
    if(any(selection(2:)>0)) valid=valid.and.depth>=0.and.stage_only
    call agree(valid)
    if(.not.valid) return
    pressure=0
    if (present(outlet_pressure)) pressure=outlet_pressure
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
    controls=[root_spacing,dt,pressure,mach,tin,pin]
    valid=all(ieee_is_finite(controls)).and.all(controls(:4)>0)
    valid=valid.and.(present(transport_rate).eqv.present(integrated_flux))
    valid=valid.and.all(pressure>=0).and.all(pressure<=1e15_dp)
    do node=1,6
      if (pressure(node)>0.and.faces(node)/='outflow') valid=.false.
    end do
    if(any(thermal_slip)) then
      if(any(thermal_slip.and.faces/='slip_wall')) valid=.false.
      if(.not.conduction.or..not.present(wall_temperature)) valid=.false.
    endif
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
    call MPI_Allreduce(controls,low,28,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high,28,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(low/=high)) return
    flags=[merge(1,0,viscosity),merge(1,0,conduction),merge(1,0,diffusion),merge(1,0,barodiffusion), &
      merge(1,0,present(wall_temperature)),merge(1,0,present(wall_velocity)), &
      merge(1,0,present(coarse_limiter)),merge(1,0,stage_only),merge(1,0,present(frame))]
    call MPI_Allreduce(flags,min_flags,size(flags),MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,max_flags,size(flags),MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(min_flags/=max_flags)) return
    if(present(frame)) then
      call MPI_Allreduce(frame%axis,axis_low,1,MPI_INTEGER,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(frame%axis,axis_high,1,MPI_INTEGER,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.axis_low/=axis_high) return
      rotation=[domain_lower,frame%angular_speed,frame%axis_point]
      call MPI_Allreduce(rotation,rotation_low,size(rotation),MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if(status/=MPI_SUCCESS) return
      call MPI_Allreduce(rotation,rotation_high,size(rotation),MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
      if(status/=MPI_SUCCESS.or.any(rotation_low/=rotation_high)) return
    endif
    options(:6)=faces;options(7)='pcm'
    if (present(coarse_limiter)) options(7)=coarse_limiter
    reference_options=options
    call MPI_Bcast(reference_options,64*7,MPI_CHARACTER,0,comm,status)
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
    if(.not.stage_only) then
      call average_down_mpi_amr_tree_3d(tree,owners,initial,nvar,comm,valid)
      if (.not.valid) return
      call recover_tree(initial,temps,valid)
      if (.not.valid) return
    end if
    call copy_local(initial,candidate,valid)
    if (.not.valid) return
    call copy_local(temps,initial_temperature,valid)
    if (.not.valid) return
    call euler_stage(candidate,temps,first_budget,first_max,first_theta,valid)
    if (.not.valid) return
    if(stage_only) then
      minimum_theta=first_theta;ok=.true.
      return
    end if
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
      call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,u,t,comm,success,internal_energy=internal_energy)
    end subroutine
    subroutine euler_stage(u,t,budget,largest,theta,success)
      type(amr_tree_state_3d), intent(inout) :: u,t
      real(dp), intent(out) :: budget(:),largest,theta
      real(dp) :: work_stamp
      logical, intent(out) :: success
      type(amr_tree_state_3d) :: packed,updated,ghosts
      type(amr_tree_flux_3d), allocatable :: flux(:)
      real(dp), allocatable :: fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:),ordinary(:,:,:,:),ordinary_ghost(:,:,:,:)
      real(dp), allocatable :: sampled_h(:,:,:),sampled_ghost_h(:,:,:)
      real(dp) :: correction(nvar),rhs(nvar),local_theta,global_theta,lower(3)
      character(len=32) :: local_faces(6)
      integer :: peer,n(3),dims(3),i,j,k,allocation_status,axis,level
      success=.false.;budget=0;largest=0;correction=0;theta=1
      allocate(packed%fields(size(tree%nodes)),flux(size(tree%nodes)),stat=allocation_status)
      success=allocation_status==0;call agree(success)
      if (.not.success) return
      do peer=1,size(tree%nodes)
        if (owners(peer)/=rank) cycle
        n=tree%nodes(peer)%shape
        allocate(packed%fields(peer)%values(np,n(1),n(2),n(3)),stat=allocation_status)
        if (allocation_status/=0) exit
        packed%fields(peer)%values(:nvar,:,:,:)=u%fields(peer)%values
        packed%fields(peer)%values(nvar+1,:,:,:)=t%fields(peer)%values(1,:,:,:)
        if(present(internal_energy)) packed%fields(peer)%values(np:np,:,:,:)=internal_energy%fields(peer)%values
      end do
      success=allocation_status==0;call agree(success)
      if (.not.success) return
      call copy_local(u,updated,success)
      if (.not.success) return
      if(stage_only) then
        do peer=1,size(tree%nodes)
          if(owners(peer)==rank) updated%fields(peer)%values=0
        end do
        if(present(coarse_start)) then
          call validate_mpi_amr_tree_local_3d(tree,owners,coarse_start,np,comm,success,maximum_depth=depth-1)
          if(.not.success) return
          do peer=1,size(tree%nodes)
            if(owners(peer)/=rank.or.tree%nodes(peer)%level>=depth) cycle
            if(.not.all(ieee_is_finite(coarse_start%fields(peer)%values))) success=.false.
            packed%fields(peer)%values=coarse_start%fields(peer)%values
          end do
          call agree(success)
          if(.not.success) return
        end if
      end if
      call fill_mpi_amr_tree_ghosts_3d(tree,owners,packed,np,faces,ghosts,comm,success, &
        viscous_wall=.true.,inflow_state=inlet,wall_velocity=wall_velocity, &
        coarse_limiter=sampling_limiter,species_count=size(species), &
        species=species,inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure, &
        depth=target_depth,coarse_future=coarse_end,coarse_fractions=coarse_fractions, &
        frame=frame,domain_lower=domain_lower,root_spacing=root_spacing, &
        stored_internal_energy=present(internal_energy),limit_internal_energy=present(internal_energy))
      if (.not.success) return
      do peer=1,size(tree%nodes)
        if (owners(peer)/=rank) cycle
        if(depth>=0.and.tree%nodes(peer)%level/=depth) cycle
        n=tree%nodes(peer)%shape;dims=n+4
        work_stamp=tree_work_start()
        associate(ghost=>ghosts%fields(peer)%values)
          ordinary=u%fields(peer)%values
          ordinary_ghost=ghost
          if(present(frame)) then
            lower=domain_lower+real(origin(:,peer),dp)*h(:,peer)
            call convert_regular_energy_3d(frame,lower,h(:,peer),ordinary,0,.false.,success)
            if(.not.success) exit
            call convert_regular_energy_3d(frame,lower,h(:,peer),ordinary_ghost,2,.false.,success)
            if(.not.success) exit
          endif
          call apply_amr_tree_pressure_ghosts_3d(tree,species,peer,faces,pressure,ordinary_ghost,success, &
            stored_internal_energy=present(internal_energy))
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
          ! Match serial MOL total-energy primitives without changing donor packing.
          if(present(internal_energy).and..not.use_total_energy) then
            sampled_h=internal_energy%fields(peer)%values(1,:,:,:)
            sampled_ghost_h=ordinary_ghost(np,:,:,:)
          endif
          local_thermal_slip=thermal_slip.and.local_faces=='slip_wall'
          call physical_transport_ghosted_fluxes_3d(species,transport,ordinary, &
            t%fields(peer)%values(1,:,:,:),ordinary_ghost(:nvar,:,:,:),ordinary_ghost(nvar+1,:,:,:), &
            n,h(:,peer),dt,local_faces, &
            viscosity,conduction,diffusion,barodiffusion,fx,fy,fz,local_theta,success,wall_temperature,wall_velocity, &
          constant_coefficients=constant_coefficients, &
            internal_energy_density=sampled_h,internal_energy_ghost=sampled_ghost_h, &
            isothermal_slip_faces=local_thermal_slip)
          theta=min(theta,local_theta)
          if (.not.success) exit
          do k=1,n(3)
            do j=1,n(2)
              do i=1,n(1)
                rhs=-(fx(:,i+1,j,k)-fx(:,i,j,k))/h(1,peer) &
                  -(fy(:,i,j+1,k)-fy(:,i,j,k))/h(2,peer) &
                  -(fz(:,i,j,k+1)-fz(:,i,j,k))/h(3,peer)
                if(stage_only) then
                  updated%fields(peer)%values(:,i,j,k)=rhs
                else
                  updated%fields(peer)%values(:,i,j,k)=u%fields(peer)%values(:,i,j,k)+dt*rhs
                end if
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
      call MPI_Allreduce(theta,global_theta,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
      if (status/=MPI_SUCCESS) then
        success=.false.;return
      end if
      theta=global_theta
      if(stage_only) then
        do peer=1,size(tree%nodes)
          if(owners(peer)/=rank) cycle
          if(.not.all(ieee_is_finite(updated%fields(peer)%values))) success=.false.
          if(depth>=0.and.tree%nodes(peer)%level/=depth) cycle
          if(.not.all(ieee_is_finite(flux(peer)%x)).or..not.all(ieee_is_finite(flux(peer)%y)).or. &
            .not.all(ieee_is_finite(flux(peer)%z))) success=.false.
        end do
        call agree(success)
        if(.not.success) return
        call move_alloc(updated%fields,transport_rate%fields)
        call move_alloc(flux,integrated_flux)
        return
      end if
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
