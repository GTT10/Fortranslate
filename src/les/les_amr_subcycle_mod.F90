module les_amr_subcycle_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  use reactive_1d_mod, only: reactive_nvar
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,amr_tree_integer_layout_3d, &
    average_down_amr_patch_tree_3d
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use les_amr_patch_tree_mod, only: advance_amr_patch_tree_les_level_3d,compute_amr_patch_tree_les_timestep_3d, &
    convert_amr_les_tree_energy_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d,reflux_amr_patch_tree_3d
  implicit none
  private
  public :: advance_amr_patch_tree_subcycle_les_3d,compute_amr_patch_tree_subcycle_les_timestep_3d
contains
  subroutine compute_amr_patch_tree_subcycle_les_timestep_3d(tree,species,state,temperature,root_spacing,cfl,faces, &
      cs,prt,sct,dt,maximum_diffusivity,ok,inflow_state,inflow_temperature,wall_temperature,wall_velocity, &
      coarse_limiter,outlet_pressure,model,ci,test_model,test_ratio,grid_ratio,harmonic_mean,inflow_mach,inflow_pressure, &
      frame,domain_lower)
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: root_spacing(3),cfl,cs,prt,sct
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(inout) :: dt,maximum_diffusivity
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: outlet_pressure(6),ci
    character(len=*), optional, intent(in) :: coarse_limiter,model
    character(len=*), optional, intent(in) :: test_model
    integer, optional, intent(in) :: test_ratio,grid_ratio
    logical, optional, intent(in) :: harmonic_mean
    call compute_amr_patch_tree_les_timestep_3d(tree,species,state,temperature,root_spacing,cfl,faces, &
      cs,prt,sct,dt,maximum_diffusivity,ok,inflow_state,inflow_temperature,wall_temperature,wall_velocity, &
      coarse_limiter,outlet_pressure,model,ci,test_model,test_ratio,grid_ratio,harmonic_mean,subcycle=.true., &
      inflow_mach=inflow_mach,inflow_pressure=inflow_pressure,frame=frame,domain_lower=domain_lower)
  end subroutine

  recursive subroutine advance_amr_patch_tree_subcycle_les_3d(tree,species,state,temperature,spacing,dt,faces, &
      cs,prt,sct,boundary_change,maximum_correction,maximum_diffusivity,ok, &
      inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter,level_steps,outlet_pressure, &
      model,ci,test_model,test_ratio,grid_ratio,harmonic_mean,inflow_mach,inflow_pressure,frame,domain_lower)
    type(rotating_frame_3d), optional, intent(in) :: frame
    real(dp), optional, intent(in) :: domain_lower(3)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: spacing(3),dt
    character(len=*), intent(in) :: faces(6)
    character(len=*), optional, intent(in) :: coarse_limiter
    real(dp), intent(in) :: cs,prt,sct
    character(len=*), optional, intent(in) :: model,test_model
    real(dp), optional, intent(in) :: ci
    integer, optional, intent(in) :: test_ratio,grid_ratio
    logical, optional, intent(in) :: harmonic_mean
    real(dp), optional, intent(in) :: wall_temperature(:),wall_velocity(:,:)
    real(dp), intent(inout) :: boundary_change(:),maximum_correction,maximum_diffusivity
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:)
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: outlet_pressure(6)
    integer(int64), optional, intent(inout) :: level_steps(:)
    type(amr_tree_state_3d) :: candidate,temps,old_endpoint,new_endpoint
    type(amr_tree_state_3d) :: ordinary,ordinary_t
    type(amr_tree_flux_3d), allocatable :: integrated(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:),counts(:),executed(:)
    integer(int64), allocatable :: rotated_steps(:)
    integer, allocatable :: ratios(:)
    real(dp), allocatable :: starts(:),finishes(:)
    real(dp) :: budget(size(boundary_change)),largest,diffusivity
    integer :: nvar,node,n(3)
    logical :: valid
    real(dp) :: rotated_budget(size(boundary_change)),rotated_max,rotated_diffusivity
    ok=.false.
    if(present(frame).neqv.present(domain_lower)) return
    if(present(frame)) then
      if(present(level_steps)) then
        if(size(level_steps)/=tree%level_count()) return
      endif
      ordinary=state;ordinary_t=temperature
      call convert_amr_les_tree_energy_3d(tree,ordinary,spacing,frame,domain_lower,.false.,valid)
      if(.not.valid) return
      allocate(rotated_steps(tree%level_count()))
      rotated_steps=0
      call advance_amr_patch_tree_subcycle_les_3d(tree,species,ordinary,ordinary_t,spacing,dt,faces, &
        cs,prt,sct,rotated_budget,rotated_max,rotated_diffusivity,valid, &
        inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter,rotated_steps,outlet_pressure, &
        model,ci,test_model,test_ratio,grid_ratio,harmonic_mean,inflow_mach,inflow_pressure)
      if(.not.valid) return
      call convert_amr_les_tree_energy_3d(tree,ordinary,spacing,frame,domain_lower,.true.,valid)
      if(.not.valid) return
      state=ordinary;temperature=ordinary_t;boundary_change=rotated_budget
      maximum_correction=rotated_max;maximum_diffusivity=rotated_diffusivity
      if(present(level_steps)) level_steps=rotated_steps
      ok=.true.;return
    endif
    call layout(tree,species,state,temperature,spacing,origin,domain,ratios,counts,valid)
    if (.not.valid) return
    nvar=reactive_nvar(size(species))
    if (size(boundary_change)/=nvar.or..not.ieee_is_finite(dt)) return
    if (dt<=0.or.dt>1e6_dp.or..not.valid_reactive_boundaries_3d(faces)) return
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
    starts=0;finishes=1;executed=0;budget=0;largest=0;diffusivity=0
    call interval(0,0._dp,1._dp,integrated,valid)
    if (.not.valid) return
    if (any(executed/=counts)) return
    state=candidate;temperature=temps;boundary_change=budget;maximum_correction=largest;maximum_diffusivity=diffusivity
    if (present(level_steps)) level_steps=executed
    ok=.true.
  contains
    recursive subroutine interval(depth,begin_time,end_time,effective,success)
      integer, intent(in) :: depth
      real(dp), intent(in) :: begin_time,end_time
      type(amr_tree_flux_3d), allocatable, intent(out) :: effective(:)
      logical, intent(out) :: success
      type(amr_tree_flux_3d), allocatable :: total(:),child_flux(:)
      real(dp) :: fraction0(depth),fraction1(depth),correction(nvar),max_change,sub_begin,sub_end,local_dt,local_diffusivity
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
      local_diffusivity=0
      call advance_amr_patch_tree_les_level_3d(tree,species,candidate,temps,depth,spacing,local_dt,faces, &
        cs,prt,sct,total,local_diffusivity,good,old_endpoint,new_endpoint,fraction0,fraction1, &
        inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter,outlet_pressure, &
        model,ci,test_model,test_ratio,grid_ratio,harmonic_mean,inflow_mach,inflow_pressure)
      if (.not.good) return
      diffusivity=max(diffusivity,local_diffusivity)
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
end module les_amr_subcycle_mod
