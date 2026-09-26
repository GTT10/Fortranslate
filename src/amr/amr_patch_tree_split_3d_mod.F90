module amr_patch_tree_split_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use amr_patch_tree_subcycle_3d_mod, only: advance_amr_patch_tree_subcycle_hydro_3d
  use amr_patch_tree_transport_subcycle_3d_mod, only: advance_amr_patch_tree_subcycle_transport_3d
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d, &
    average_down_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  use amr_patch_tree_chemistry_3d_mod, only: advance_amr_patch_tree_chemistry_3d,react_sdc_amr_patch_tree_3d
  use amr_patch_tree_hydro_3d_mod, only: advance_amr_patch_tree_hydro_3d,compute_sdc_amr_patch_tree_hydro_3d
  use amr_patch_tree_transport_3d_mod, only: advance_amr_patch_tree_transport_3d, &
    compute_sdc_amr_patch_tree_transport_3d
  use reactive_3d_mod, only: compute_reactive_external_source_3d,recover_reactive_temperatures_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d,reflux_amr_patch_tree_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  implicit none
  private
  public :: advance_amr_patch_tree_split_3d,advance_amr_patch_tree_sdc_3d
contains
  ! Whole-tree or level-local SDC. Level mode exports fluxes without reflux;
  ! a recursive caller owns child intervals and final synchronization.
  subroutine advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,state,temperature,spacing,dt,faces, &
      iterations,chemistry,transport_flags,force,rtol,atol,boundary_change,source_change,maximum_correction,minimum_theta, &
      ok,reaction_source,inflow_state,inflow_temperature,wall_temperature,wall_velocity,outlet_pressure, &
      inflow_mach,inflow_pressure,difmag,target_depth,coarse_start,coarse_end,fractions_start,fractions_end, &
      source_start,source_end,integrated_flux,old_predictor,new_predictor)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: spacing(3),dt,force(3),rtol,atol
    character(len=*), intent(in) :: faces(6)
    integer, intent(in) :: iterations
    logical, intent(in) :: chemistry,transport_flags(4)
    real(dp), intent(inout) :: boundary_change(:),source_change(:),maximum_correction,minimum_theta
    logical, intent(out) :: ok
    type(amr_tree_state_3d), optional, intent(inout) :: reaction_source
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6),difmag
    integer, optional, intent(in) :: target_depth
    type(amr_tree_state_3d), optional, intent(in) :: coarse_start,coarse_end,source_start,source_end
    real(dp), optional, intent(in) :: fractions_start(:),fractions_end(:)
    type(amr_tree_flux_3d), allocatable, optional, intent(inout) :: integrated_flux(:)
    type(amr_tree_state_3d), optional, intent(inout) :: old_predictor,new_predictor
    type(amr_tree_state_3d) :: old,ot,candidate,ct,old_ext,new_ext,old_tr,new_tr,hydro,rates,ir,op,np
    real(dp), allocatable :: midpoint(:)
    type(amr_tree_flux_3d), allocatable :: flux(:),old_flux(:),new_flux(:)
    real(dp) :: budget(size(boundary_change)),volume(size(source_change)),correction(size(boundary_change))
    real(dp) :: first(size(source_change)),last(size(source_change)),largest,theta,stage_theta
    integer :: nv,node,n(3),iteration,c,depth
    logical :: valid,do_transport
    ok=.false.;nv=5+size(species)
    if(.not.state%is_valid(tree).or..not.temperature%is_valid(tree)) return
    if(size(state%fields(1)%values,1)/=nv.or.size(temperature%fields(1)%values,1)/=1) return
    if(size(boundary_change)/=nv.or.size(source_change)/=nv.or.iterations<1) return
    if(.not.all(ieee_is_finite([spacing,dt,force,rtol,atol]))) return
    if(minval(spacing)<=0.or.min(dt,rtol,atol)<=0.or..not.valid_reactive_boundaries_3d(faces)) return
    if(transport_flags(4).and..not.transport_flags(3)) return
    depth=-1
    if(present(target_depth)) then
      depth=target_depth
      if(depth<0.or.depth>=tree%level_count()) return
      if(.not.present(integrated_flux)) return
      if(depth>0.and..not.present(coarse_start)) return
    end if
    if(present(coarse_start).neqv.present(coarse_end)) return
    if(present(coarse_end).neqv.present(fractions_start)) return
    if(present(fractions_start).neqv.present(fractions_end)) return
    if(present(coarse_start).neqv.present(source_start)) return
    if(present(source_start).neqv.present(source_end)) return
    if(present(coarse_start)) then
      if(depth<0) return
      if(size(fractions_start)/=depth.or.size(fractions_end)/=depth) return
      if(.not.all(ieee_is_finite(fractions_start)).or..not.all(ieee_is_finite(fractions_end))) return
      if(any(fractions_start<0).or.any(fractions_end>1).or.any(fractions_end<fractions_start)) return
      midpoint=.5_dp*(fractions_start+fractions_end)
    end if
    do_transport=any(transport_flags(:3));theta=1
    old=state;ot=temperature
    if(depth<0) then
      call average_down_amr_patch_tree_3d(tree,old,valid)
      if(.not.valid) return
    end if
    call recover(old,ot,valid)
    if(.not.valid) return
    allocate(ir%fields(size(tree%nodes)))
    do node=1,size(tree%nodes)
      n=tree%nodes(node)%shape
      allocate(ir%fields(node)%values(nv-4,n(1),n(2),n(3)));ir%fields(node)%values=0
    end do
    if(present(reaction_source)) then
      if(.not.chemistry.or..not.reaction_source%is_valid(tree)) return
      if(size(reaction_source%fields(1)%values,1)/=nv-4) return
      ir=reaction_source
    end if
    call sources(old,ot,old_ext,old_tr,old_flux,valid,fractions_start)
    if(.not.valid) return
    call predictor(old_ext,old_tr,op)
    new_ext=old_ext;new_tr=old_tr;candidate=old;ct=ot
    do iteration=1,iterations
      call average_rates()
      do node=1,size(tree%nodes)
        if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
        rates%fields(node)%values(6:,:,:,:)=rates%fields(node)%values(6:,:,:,:)+ir%fields(node)%values(:nv-5,:,:,:)
        rates%fields(node)%values(5,:,:,:)=rates%fields(node)%values(5,:,:,:)+ir%fields(node)%values(nv-4,:,:,:)
      end do
      call compute_sdc_amr_patch_tree_hydro_3d(tree,species,old,ot,rates,spacing,dt,faces,hydro,flux,valid, &
        inflow_state,inflow_temperature,outlet_pressure,inflow_mach,inflow_pressure,difmag, &
        target_depth,coarse_start,coarse_end,fractions_start,source_start,source_end,midpoint)
      if(.not.valid) return
      do node=1,size(tree%nodes)
        if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
        candidate%fields(node)%values=old%fields(node)%values+dt*(rates%fields(node)%values+hydro%fields(node)%values)
      end do
      call recover(candidate,ct,valid)
      if(.not.valid) return
      call sources(candidate,ct,new_ext,new_tr,new_flux,valid,fractions_end)
      if(.not.valid) return
      call average_rates()
      do node=1,size(tree%nodes)
        if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
        rates%fields(node)%values=rates%fields(node)%values+hydro%fields(node)%values
      end do
      if(chemistry) then
        call react_sdc_amr_patch_tree_3d(tree,species,reactions,old,ot,rates,dt,rtol,atol,candidate,ct,ir,valid,target_depth)
        if(.not.valid) return
      else
        do node=1,size(tree%nodes)
        if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
          candidate%fields(node)%values=old%fields(node)%values+dt*rates%fields(node)%values
        end do
        call recover(candidate,ct,valid)
        if(.not.valid) return
      end if
    end do
    if(do_transport) then
      do node=1,size(tree%nodes)
        if(depth>=0.and.tree%nodes(node)%level/=depth) cycle
        flux(node)%x=flux(node)%x+.5_dp*(old_flux(node)%x+new_flux(node)%x)
        flux(node)%y=flux(node)%y+.5_dp*(old_flux(node)%y+new_flux(node)%y)
        flux(node)%z=flux(node)%z+.5_dp*(old_flux(node)%z+new_flux(node)%z)
      end do
    end if
    budget=0;largest=0
    if(depth<=0) then
      n=tree%nodes(1)%shape
      do c=1,nv
        budget(c)=(sum(flux(1)%x(c,1,:,:))-sum(flux(1)%x(c,n(1)+1,:,:)))*spacing(2)*spacing(3)+ &
          (sum(flux(1)%y(c,:,1,:))-sum(flux(1)%y(c,:,n(2)+1,:)))*spacing(1)*spacing(3)+ &
          (sum(flux(1)%z(c,:,:,1))-sum(flux(1)%z(c,:,:,n(3)+1)))*spacing(1)*spacing(2)
      end do
    end if
    if(depth<0) then
      call reflux_amr_patch_tree_3d(tree,candidate,flux,spacing,faces,correction,largest,valid)
      if(.not.valid) return
      budget=budget+correction
      call recover(candidate,ct,valid)
      if(.not.valid) return
    end if
    call composite_integrals_amr_patch_tree_3d(tree,old_ext,spacing,first,valid)
    if(.not.valid) return
    call composite_integrals_amr_patch_tree_3d(tree,new_ext,spacing,last,valid)
    if(.not.valid) return
    volume=.5_dp*dt*(first+last)
    if(.not.all(ieee_is_finite(budget)).or..not.all(ieee_is_finite(volume))) return
    if(present(old_predictor)) then
      if(.not.op%is_valid(tree)) return
    end if
    if(present(new_predictor)) then
      call predictor(new_ext,new_tr,np)
      if(.not.np%is_valid(tree)) return
    end if
    state=candidate;temperature=ct;boundary_change=budget;source_change=volume
    maximum_correction=largest;minimum_theta=theta
    if(present(reaction_source)) reaction_source=ir
    if(present(integrated_flux)) integrated_flux=flux
    if(present(old_predictor)) old_predictor=op
    if(present(new_predictor)) new_predictor=np
    ok=.true.
  contains
    subroutine predictor(external,diffusive,result)
      type(amr_tree_state_3d), intent(in) :: external,diffusive
      type(amr_tree_state_3d), intent(out) :: result
      integer :: peer
      result=external
      do peer=1,size(tree%nodes)
        if(depth>=0.and.tree%nodes(peer)%level/=depth) cycle
        result%fields(peer)%values=result%fields(peer)%values+diffusive%fields(peer)%values
        result%fields(peer)%values(6:,:,:,:)=result%fields(peer)%values(6:,:,:,:)+ir%fields(peer)%values(:nv-5,:,:,:)
        result%fields(peer)%values(5,:,:,:)=result%fields(peer)%values(5,:,:,:)+ir%fields(peer)%values(nv-4,:,:,:)
      end do
    end subroutine
    subroutine average_rates()
      integer :: peer
      rates=old_ext
      do peer=1,size(tree%nodes)
        rates%fields(peer)%values=.5_dp*new_ext%fields(peer)%values+.5_dp*old_ext%fields(peer)%values
        rates%fields(peer)%values=rates%fields(peer)%values+.5_dp*new_tr%fields(peer)%values
        rates%fields(peer)%values=rates%fields(peer)%values+.5_dp*old_tr%fields(peer)%values
      end do
    end subroutine
    subroutine recover(u,t,good)
      type(amr_tree_state_3d), intent(in) :: u
      type(amr_tree_state_3d), intent(inout) :: t
      logical, intent(out) :: good
      real(dp), allocatable :: recovered(:,:,:)
      integer :: peer,dims(3)
      do peer=1,size(tree%nodes)
        if(depth>=0.and.tree%nodes(peer)%level/=depth) cycle
        dims=tree%nodes(peer)%shape
        allocate(recovered(dims(1),dims(2),dims(3)))
        call recover_reactive_temperatures_3d(species,u%fields(peer)%values,t%fields(peer)%values(1,:,:,:), &
          dims(1),dims(2),dims(3),recovered,good)
        if(.not.good) return
        t%fields(peer)%values(1,:,:,:)=recovered
        deallocate(recovered)
      end do
      good=.true.
    end subroutine
    subroutine sources(u,t,external,diffusive,face_flux,good,fractions)
      type(amr_tree_state_3d), intent(in) :: u,t
      type(amr_tree_state_3d), intent(inout) :: external,diffusive
      type(amr_tree_flux_3d), allocatable, intent(inout) :: face_flux(:)
      logical, intent(out) :: good
      real(dp), optional, intent(in) :: fractions(:)
      integer :: peer
      external=u;diffusive=u
      do peer=1,size(tree%nodes)
        external%fields(peer)%values=0;diffusive%fields(peer)%values=0
        if(depth>=0.and.tree%nodes(peer)%level/=depth) cycle
        call compute_reactive_external_source_3d(u%fields(peer)%values,force,external%fields(peer)%values,good)
        if(.not.good) return
        diffusive%fields(peer)%values=0
      end do
      if(do_transport) then
        stage_theta=1
        call compute_sdc_amr_patch_tree_transport_3d(tree,species,transport,u,t,spacing,dt,faces,transport_flags, &
          diffusive,face_flux,stage_theta,good,inflow_state,inflow_temperature,wall_temperature,wall_velocity, &
          outlet_pressure,inflow_mach,inflow_pressure,target_depth,coarse_start,coarse_end,fractions)
        if(.not.good) return
        theta=min(theta,stage_theta)
      end if
      good=.true.
    end subroutine
  end subroutine advance_amr_patch_tree_sdc_3d

  ! R/2 - T/2 - H - T/2 - R/2, over one coarse interval, with optional recursive flow/transport steps.
  ! Disabled operators are skipped. All stage outputs remain private until
  ! the complete composition succeeds, including the final chemistry solve.
  subroutine advance_amr_patch_tree_split_3d(tree,species,reactions,transport,state,temperature,spacing,dt,faces, &
      reconstruction,limiter,solver,chemistry,viscosity,conduction,diffusion,barodiffusion,rtol,atol, &
      boundary_change,maximum_correction,minimum_theta,ok,chemistry_integrator, &
      inflow_state,inflow_temperature,wall_temperature,wall_velocity,subcycling,outlet_pressure,inflow_mach,inflow_pressure, &
      difmag)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: spacing(3),dt,rtol,atol
    character(len=*), intent(in) :: faces(6),reconstruction,limiter,solver
    logical, intent(in) :: chemistry,viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(inout) :: boundary_change(:),maximum_correction,minimum_theta
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: chemistry_integrator
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    real(dp), optional, intent(in) :: outlet_pressure(6),inflow_mach(6),inflow_pressure(6)
    real(dp), optional, intent(in) :: difmag
    character(len=:), allocatable :: transport_limiter
    type(amr_tree_state_3d) :: candidate,temps
    real(dp) :: budget(size(boundary_change)),stage_budget(size(boundary_change)),largest,stage_max,theta,stage_theta
    logical, optional, intent(in) :: subcycling
    logical :: valid,transport_enabled,recursive_steps
    ok=.false.
    if (.not.ieee_is_finite(dt)) return
    if (dt<=0.or.(barodiffusion.and..not.diffusion)) return
    if (reconstruction=='characteristic_plm'.or.reconstruction=='characteristic_ppm') transport_limiter=limiter
    recursive_steps=.false.
    if (present(subcycling)) recursive_steps=subcycling
    candidate=state;temps=temperature
    budget=0;largest=0;theta=1
    transport_enabled=viscosity.or.conduction.or.diffusion
    if (chemistry) then
      call reaction_half(valid)
      if (.not.valid) return
    end if
    if (transport_enabled) then
      call transport_half(valid)
      if (.not.valid) return
    end if
    stage_budget=0;stage_max=0
    if (recursive_steps) then
      call advance_amr_patch_tree_subcycle_hydro_3d(tree,species,candidate,temps,spacing,dt,faces, &
        reconstruction,limiter,solver,stage_budget,stage_max,valid,inflow_state,inflow_temperature, &
        outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure,difmag=difmag)
    else
      call advance_amr_patch_tree_hydro_3d(tree,species,candidate,temps,spacing,dt,faces,reconstruction,limiter,solver, &
        stage_budget,stage_max,valid,inflow_state,inflow_temperature, &
        outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure,difmag=difmag)
    end if
    if (.not.valid) return
    budget=budget+stage_budget;largest=max(largest,stage_max)
    if (transport_enabled) then
      call transport_half(valid)
      if (.not.valid) return
    end if
    if (chemistry) then
      call reaction_half(valid)
      if (.not.valid) return
    end if
    state=candidate;temperature=temps
    boundary_change=budget;maximum_correction=largest;minimum_theta=theta;ok=.true.
  contains
    subroutine reaction_half(success)
      logical, intent(out) :: success
      call advance_amr_patch_tree_chemistry_3d(tree,species,reactions,candidate,temps,.5_dp*dt,rtol,atol, &
        success,chemistry_integrator)
    end subroutine
    subroutine transport_half(success)
      logical, intent(out) :: success
      stage_budget=0;stage_max=0;stage_theta=1
      if (recursive_steps) then
        call advance_amr_patch_tree_subcycle_transport_3d(tree,species,transport,candidate,temps,spacing,.5_dp*dt,faces, &
          viscosity,conduction,diffusion,barodiffusion,stage_budget,stage_max,stage_theta,success, &
          inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter=transport_limiter, &
          outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
      else
        call advance_amr_patch_tree_transport_3d(tree,species,transport,candidate,temps,spacing,.5_dp*dt,faces, &
          viscosity,conduction,diffusion,barodiffusion,stage_budget,stage_max,stage_theta,success, &
          inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter=transport_limiter, &
          outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
      end if
      if (.not.success) return
      budget=budget+stage_budget;largest=max(largest,stage_max);theta=min(theta,stage_theta)
    end subroutine
  end subroutine advance_amr_patch_tree_split_3d
end module amr_patch_tree_split_3d_mod
