module amr_patch_tree_transport_subcycle_3d_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use reactive_transport_3d_mod, only: reactive_transport_timestep_3d,reactive_isothermal_timestep_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  use reactive_1d_mod, only: reactive_nvar
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d, &
    average_down_amr_patch_tree_3d
  use amr_patch_tree_transport_3d_mod, only: advance_amr_patch_tree_transport_level_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d,reflux_amr_patch_tree_3d
  implicit none
  private
  public :: compute_amr_patch_tree_subcycle_transport_timestep_3d,advance_amr_patch_tree_subcycle_transport_3d
contains
  subroutine compute_amr_patch_tree_subcycle_transport_timestep_3d(tree,species,transport,state,temperature, &
      root_spacing,cfl,viscosity,conduction,diffusion,dt,maximum_diffusivity,ok,faces,wall_temperature, &
      constant_coefficients,internal_energy,isothermal_slip_faces)
    type(amr_tree_state_3d), optional, intent(in) :: internal_energy
    real(dp), allocatable :: sampled_h(:,:,:)
    real(dp), optional, intent(in) :: constant_coefficients(:)
    type(amr_patch_tree_3d), intent(in) :: tree
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
    integer :: node,n(3),axis,level
    logical, optional, intent(in) :: isothermal_slip_faces(6)
    logical :: thermal_slip(6)
    logical :: valid
    ok=.false.
    thermal_slip=.false.
    if(present(isothermal_slip_faces)) thermal_slip=isothermal_slip_faces
    if(any(thermal_slip)) then
      if(.not.present(faces).or..not.present(wall_temperature).or..not.conduction) return
      if(size(wall_temperature)/=6) return
      if(any(thermal_slip.and.faces/='slip_wall')) return
      if(any(thermal_slip.and.(.not.ieee_is_finite(wall_temperature).or.wall_temperature<=0))) return
    endif
    call layout(tree,species,state,temperature,root_spacing,origin,domain,ratios,counts,valid)
    if (.not.valid) return
    if (present(faces)) then
      if (.not.valid_reactive_boundaries_3d(faces)) return
      if (any(faces=='isothermal_wall').or.any(thermal_slip)) then
        if (.not.conduction.or..not.present(wall_temperature)) return
        if (size(wall_temperature)/=6) return
      end if
    end if
    allocate(h(3,size(tree%nodes)))
    do node=1,size(tree%nodes)
      h(:,node)=root_spacing*real(domain(:,1),dp)/real(domain(:,tree%nodes(node)%level+1),dp)
    end do
    if(present(internal_energy)) then
      if(.not.internal_energy%is_valid(tree)) return
      do node=1,size(tree%nodes)
        if(size(internal_energy%fields(node)%values,1)/=1) return
        if(.not.all(ieee_is_finite(internal_energy%fields(node)%values))) return
      enddo
    endif
    candidate=huge(1.0_dp);diffusivity=0
    do node=1,size(tree%nodes)
      n=tree%nodes(node)%shape;level=tree%nodes(node)%level+1
      if(present(internal_energy)) sampled_h=internal_energy%fields(node)%values(1,:,:,:)
      call reactive_transport_timestep_3d(species,transport,state%fields(node)%values,temperature%fields(node)%values(1,:,:,:), &
        n(1),n(2),n(3),h(1,node),h(2,node),h(3,node),cfl,viscosity,conduction,diffusion,local_dt,local_diffusivity,valid, &
          constant_coefficients=constant_coefficients,internal_energy_density=sampled_h)
      if (.not.valid) return
      if (local_dt<=huge(1._dp)/real(counts(level),dp)) candidate=min(candidate,local_dt*real(counts(level),dp))
      diffusivity=max(diffusivity,local_diffusivity)
      if (present(faces)) then
        if (any(faces=='isothermal_wall').or.any(thermal_slip)) then
          local_faces='outflow';level=tree%nodes(node)%level+1
          do axis=1,3
            if (origin(axis,node)==0.and.(faces(2*axis-1)=='isothermal_wall'.or.thermal_slip(2*axis-1))) &
              local_faces(2*axis-1)=faces(2*axis-1)
            if (origin(axis,node)+int(n(axis),int64)==domain(axis,level).and. &
                (faces(2*axis)=='isothermal_wall'.or.thermal_slip(2*axis))) &
              local_faces(2*axis)=faces(2*axis)
          end do
          call reactive_isothermal_timestep_3d(species,transport,state%fields(node)%values, &
            temperature%fields(node)%values(1,:,:,:),h(:,node),local_faces,wall_temperature,cfl,wall_dt,valid, &
          constant_coefficients=constant_coefficients,internal_energy_density=sampled_h, &
          isothermal_slip_faces=thermal_slip.and.local_faces=='slip_wall')
          if (.not.valid) return
          if (wall_dt<=huge(1._dp)/real(counts(level),dp)) candidate=min(candidate,wall_dt*real(counts(level),dp))
        end if
      end if
    end do
    dt=candidate;maximum_diffusivity=diffusivity;ok=.true.
  end subroutine

  subroutine advance_amr_patch_tree_subcycle_transport_3d(tree,species,transport,state,temperature,spacing,dt,faces, &
      viscosity,conduction,diffusion,barodiffusion,boundary_change,maximum_correction,minimum_theta,ok, &
      inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter,level_steps,outlet_pressure, &
      inflow_mach,inflow_pressure)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: spacing(3),dt
    character(len=*), intent(in) :: faces(6)
    character(len=*), optional, intent(in) :: coarse_limiter
    logical, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), optional, intent(in) :: wall_temperature(:),wall_velocity(:,:)
    real(dp), intent(inout) :: boundary_change(:),maximum_correction,minimum_theta
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:)
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6)
    integer(int64), optional, intent(inout) :: level_steps(:)
    type(amr_tree_state_3d) :: candidate,temps,old_endpoint,new_endpoint
    type(amr_tree_flux_3d), allocatable :: integrated(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:),counts(:),executed(:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: starts(:),finishes(:)
    real(dp) :: budget(size(boundary_change)),largest,theta
    integer :: nvar,node,n(3)
    logical :: valid
    ok=.false.
    call layout(tree,species,state,temperature,spacing,origin,domain,ratios,counts,valid)
    if (.not.valid) return
    nvar=reactive_nvar(size(species))
    if (size(boundary_change)/=nvar.or..not.ieee_is_finite(dt)) return
    if (dt<=0) return
    if (present(level_steps)) then
      if (size(level_steps)/=size(counts)) return
    end if
    candidate=state;temps=temperature
    call average_down_amr_patch_tree_3d(tree,candidate,valid)
    if (.not.valid) return
    call recover_from(0,valid)
    if (.not.valid) return
    allocate(old_endpoint%fields(size(tree%nodes)),new_endpoint%fields(size(tree%nodes)))
    do node=1,size(tree%nodes)
      n=tree%nodes(node)%shape
      allocate(old_endpoint%fields(node)%values(nvar+1,n(1),n(2),n(3)))
      allocate(new_endpoint%fields(node)%values(nvar+1,n(1),n(2),n(3)))
    end do
    allocate(starts(size(counts)),finishes(size(counts)),executed(size(counts)))
    starts=0;finishes=1;executed=0;budget=0;largest=0;theta=1
    call interval(0,0._dp,1._dp,integrated,valid)
    if (.not.valid) return
    if (any(executed/=counts)) return
    state=candidate;temperature=temps;boundary_change=budget;maximum_correction=largest;minimum_theta=theta
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
      integer :: peer,ancestor
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
      call advance_amr_patch_tree_transport_level_3d(tree,species,transport,candidate,temps,depth,spacing,local_dt,faces, &
        viscosity,conduction,diffusion,barodiffusion,total,local_theta,good,old_endpoint,new_endpoint,fraction0,fraction1, &
        inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter,outlet_pressure,inflow_mach,inflow_pressure)
      if (.not.good) return
      theta=min(theta,local_theta)
      executed(depth+1)=executed(depth+1)+1
      call pack_depth(depth,new_endpoint)
      if (depth==0) call root_budget(total(1))
      if (depth+1<size(counts)) then
        ratio=int(ratios(depth+2),int64)**2
        sub_begin=begin_time
        do substep=1,ratio
          sub_end=begin_time+(end_time-begin_time)*real(substep,dp)/real(ratio,dp)
          if (substep==ratio) sub_end=end_time
          call interval(depth+1,sub_begin,sub_end,child_flux,good)
          if (.not.good) return
          do peer=1,size(tree%nodes)
            if (tree%nodes(peer)%level<=depth) cycle
            if (.not.allocated(total(peer)%x)) then
              total(peer)=child_flux(peer)
            else
              total(peer)%x=total(peer)%x+child_flux(peer)%x
              total(peer)%y=total(peer)%y+child_flux(peer)%y
              total(peer)%z=total(peer)%z+child_flux(peer)%z
            end if
          end do
          sub_begin=sub_end
        end do
      end if
      correction=0;max_change=0
      call reflux_amr_patch_tree_3d(tree,candidate,total,spacing,faces,correction,max_change,good, &
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
        if (tree%nodes(peer)%level/=depth) cycle
        packed%fields(peer)%values(:nvar,:,:,:)=candidate%fields(peer)%values
        packed%fields(peer)%values(nvar+1,:,:,:)=temps%fields(peer)%values(1,:,:,:)
      end do
    end subroutine
    subroutine recover_from(depth,success)
      integer, intent(in) :: depth
      logical, intent(out) :: success
      real(dp), allocatable :: recovered(:,:,:)
      integer :: peer,dims(3)
      success=.false.
      do peer=1,size(tree%nodes)
        if (tree%nodes(peer)%level<depth) cycle
        dims=tree%nodes(peer)%shape;allocate(recovered(dims(1),dims(2),dims(3)))
        call recover_reactive_temperatures_3d(species,candidate%fields(peer)%values,temps%fields(peer)%values(1,:,:,:), &
          dims(1),dims(2),dims(3),recovered,success)
        if (.not.success) return
        temps%fields(peer)%values(1,:,:,:)=recovered;deallocate(recovered)
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

  subroutine layout(tree,species,state,temperature,spacing,origin,domain,ratios,counts,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: spacing(3)
    integer(int64), allocatable, intent(out) :: origin(:,:),domain(:,:),counts(:)
    integer, allocatable, intent(out) :: ratios(:)
    logical, intent(out) :: ok
    integer :: level
    integer(int64) :: factor
    ok=.false.
    if (.not.state%is_valid(tree).or..not.temperature%is_valid(tree)) return
    if (size(state%fields(1)%values,1)/=reactive_nvar(size(species))) return
    if (size(temperature%fields(1)%values,1)/=1) return
    if (.not.all(ieee_is_finite(spacing))) return
    if (any(spacing<=0)) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
    if (.not.ok) return
    ok=.false.;allocate(counts(size(ratios)));counts(1)=1
    do level=2,size(ratios)
      factor=int(ratios(level),int64)**2
      if (counts(level-1)>huge(1_int64)/factor) return
      counts(level)=counts(level-1)*factor
    end do
    ok=.true.
  end subroutine
end module amr_patch_tree_transport_subcycle_3d_mod
