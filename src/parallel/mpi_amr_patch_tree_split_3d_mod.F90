module mpi_amr_patch_tree_split_3d_mod
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use mpi_amr_patch_tree_chemistry_3d_mod, only: advance_mpi_amr_tree_chemistry_3d,clone_mpi_amr_tree_inputs_3d
  use mpi_amr_patch_tree_hydro_3d_mod, only: advance_mpi_amr_tree_hydro_3d
  use mpi_amr_patch_tree_transport_3d_mod, only: advance_mpi_amr_tree_transport_3d
  use mpi_amr_patch_tree_subcycle_3d_mod, only: advance_mpi_amr_tree_subcycle_hydro_3d
  use mpi_amr_patch_tree_transport_subcycle_3d_mod, only: advance_mpi_amr_tree_subcycle_transport_3d
  use reactive_3d_mod, only: compute_reactive_external_source_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use amr_patch_tree_remap_3d_mod, only: covered
  use mpi_amr_patch_tree_cells_3d_mod, only: validate_mpi_amr_tree_local_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: recover_mpi_amr_tree_temperatures_3d
  use mpi_amr_patch_tree_reflux_3d_mod, only: reflux_mpi_amr_tree_3d
  use mpi_amr_patch_tree_hydro_3d_mod, only: compute_sdc_mpi_amr_tree_hydro_3d
  use mpi_amr_patch_tree_transport_3d_mod, only: compute_sdc_mpi_amr_tree_transport_3d
  use mpi_amr_patch_tree_chemistry_3d_mod, only: react_sdc_mpi_amr_tree_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  implicit none
  private
  public :: advance_mpi_amr_tree_split_3d,advance_mpi_amr_tree_sdc_3d
contains
  ! Common-step SDC: every patch finishes a stage before the next stage starts.
  subroutine advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,state,temperature,spacing,dt,faces, &
      iterations,chemistry,transport_flags,force,rtol,atol,boundary_change,source_change,maximum_correction,minimum_theta, &
      comm,ok,reaction_source,inflow_state,inflow_temperature,wall_temperature,wall_velocity,outlet_pressure, &
      inflow_mach,inflow_pressure,difmag)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
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
    type(amr_tree_state_3d) :: old,ot,candidate,ct,old_ext,new_ext,old_tr,new_tr,hydro,rates,ir
    type(amr_tree_flux_3d), allocatable :: flux(:),old_flux(:),new_flux(:)
    real(dp) :: budget(size(boundary_change)),volume(size(source_change)),correction(size(boundary_change))
    real(dp) :: first(size(source_change)),last(size(source_change)),largest,theta,stage_theta
    integer :: nv,node,n(3),iteration,c,rank,status,ios
    integer :: switches(9),low_switches(9),high_switches(9)
    real(dp) :: controls(9),low_controls(9),high_controls(9)
    logical :: valid,do_transport
    ok=.false.;nv=5+size(species)
    controls=[spacing,dt,force,rtol,atol]
    valid=size(boundary_change)==nv.and.size(source_change)==nv.and.iterations>0
    valid=valid.and.all(ieee_is_finite(controls))
    valid=valid.and.minval(spacing)>0.and.min(dt,rtol,atol)>0.and.valid_reactive_boundaries_3d(faces)
    valid=valid.and..not.(transport_flags(4).and..not.transport_flags(3))
    valid=valid.and.(chemistry.or..not.present(reaction_source))
    call agree(valid)
    if(.not.valid) return
    switches=[iterations,merge(1,0,chemistry),merge(1,0,transport_flags),merge(1,0,present(reaction_source)), &
      size(species),size(reactions)]
    call MPI_Allreduce(switches,low_switches,9,MPI_INTEGER,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(switches,high_switches,9,MPI_INTEGER,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(low_switches/=high_switches)) return
    call MPI_Allreduce(controls,low_controls,9,MPI_DOUBLE_PRECISION,MPI_MIN,comm,status)
    if(status/=MPI_SUCCESS) return
    call MPI_Allreduce(controls,high_controls,9,MPI_DOUBLE_PRECISION,MPI_MAX,comm,status)
    if(status/=MPI_SUCCESS) return
    if(any(low_controls/=high_controls)) return
    do_transport=any(transport_flags(:3));theta=1
    call clone_mpi_amr_tree_inputs_3d(tree,owners,size(species),state,temperature,old,ot,comm,rank,valid)
    if(.not.valid) return
    call average_down_mpi_amr_tree_3d(tree,owners,old,nv,comm,valid)
    if(.not.valid) return
    call recover(old,ot,valid)
    if(.not.valid) return
    allocate(ir%fields(size(tree%nodes)),stat=ios)
    if(ios==0) then
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        n=tree%nodes(node)%shape
        allocate(ir%fields(node)%values(nv-4,n(1),n(2),n(3)),stat=ios)
        if(ios/=0) exit
        ir%fields(node)%values=0
      end do
    end if
    valid=ios==0;call agree(valid)
    if(.not.valid) return
    if(present(reaction_source)) then
      call validate_mpi_amr_tree_local_3d(tree,owners,reaction_source,nv-4,comm,valid)
      if(.not.valid) return
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        if(.not.all(ieee_is_finite(reaction_source%fields(node)%values))) valid=.false.
      end do
      call agree(valid)
      if(.not.valid) return
      ir=reaction_source
    end if
    call sources(old,ot,old_ext,old_tr,old_flux,valid)
    if(.not.valid) return
    new_ext=old_ext;new_tr=old_tr;candidate=old;ct=ot
    do iteration=1,iterations
      call average_rates()
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        rates%fields(node)%values(6:,:,:,:)=rates%fields(node)%values(6:,:,:,:)+ir%fields(node)%values(:nv-5,:,:,:)
        rates%fields(node)%values(5,:,:,:)=rates%fields(node)%values(5,:,:,:)+ir%fields(node)%values(nv-4,:,:,:)
      end do
      call compute_sdc_mpi_amr_tree_hydro_3d(tree,owners,species,old,ot,rates,spacing,dt,faces,hydro,flux,comm,valid, &
        inflow_state,inflow_temperature,outlet_pressure,inflow_mach,inflow_pressure,difmag)
      if(.not.valid) return
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        candidate%fields(node)%values=old%fields(node)%values+dt*(rates%fields(node)%values+hydro%fields(node)%values)
      end do
      call recover(candidate,ct,valid)
      if(.not.valid) return
      call sources(candidate,ct,new_ext,new_tr,new_flux,valid)
      if(.not.valid) return
      call average_rates()
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        rates%fields(node)%values=rates%fields(node)%values+hydro%fields(node)%values
      end do
      if(chemistry) then
        call react_sdc_mpi_amr_tree_3d(tree,owners,species,reactions,old,ot,rates,dt,rtol,atol,candidate,ct,ir,comm,valid)
        if(.not.valid) return
      else
        do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
          candidate%fields(node)%values=old%fields(node)%values+dt*rates%fields(node)%values
        end do
        call recover(candidate,ct,valid)
        if(.not.valid) return
      end if
    end do
    if(do_transport) then
      do node=1,size(tree%nodes)
        if(owners(node)/=rank) cycle
        flux(node)%x=flux(node)%x+.5_dp*(old_flux(node)%x+new_flux(node)%x)
        flux(node)%y=flux(node)%y+.5_dp*(old_flux(node)%y+new_flux(node)%y)
        flux(node)%z=flux(node)%z+.5_dp*(old_flux(node)%z+new_flux(node)%z)
      end do
    end if
    budget=0
    if(rank==owners(1)) then
    n=tree%nodes(1)%shape
    do c=1,nv
      budget(c)=(sum(flux(1)%x(c,1,:,:))-sum(flux(1)%x(c,n(1)+1,:,:)))*spacing(2)*spacing(3)+ &
        (sum(flux(1)%y(c,:,1,:))-sum(flux(1)%y(c,:,n(2)+1,:)))*spacing(1)*spacing(3)+ &
        (sum(flux(1)%z(c,:,:,1))-sum(flux(1)%z(c,:,:,n(3)+1)))*spacing(1)*spacing(2)
    end do
    end if
    call MPI_Bcast(budget,nv,MPI_DOUBLE_PRECISION,owners(1),comm,status)
    if(status/=MPI_SUCCESS) return
    call reflux_mpi_amr_tree_3d(tree,owners,candidate,nv,flux,spacing,faces,correction,largest,comm,valid)
    if(.not.valid) return
    budget=budget+correction
    call recover(candidate,ct,valid)
    if(.not.valid) return
    call inventory(old_ext,first,valid)
    if(.not.valid) return
    call inventory(new_ext,last,valid)
    if(.not.valid) return
    volume=.5_dp*dt*(first+last)
    if(.not.all(ieee_is_finite(budget)).or..not.all(ieee_is_finite(volume))) return
    state=candidate;temperature=ct;boundary_change=budget;source_change=volume
    maximum_correction=largest;minimum_theta=theta
    if(present(reaction_source)) reaction_source=ir
    ok=.true.
  contains
    subroutine average_rates()
      integer :: peer
      rates=old_ext
      do peer=1,size(tree%nodes)
        if(owners(peer)/=rank) cycle
        rates%fields(peer)%values=.5_dp*new_ext%fields(peer)%values+.5_dp*old_ext%fields(peer)%values
        rates%fields(peer)%values=rates%fields(peer)%values+.5_dp*new_tr%fields(peer)%values
        rates%fields(peer)%values=rates%fields(peer)%values+.5_dp*old_tr%fields(peer)%values
      end do
    end subroutine
    subroutine recover(u,t,good)
      type(amr_tree_state_3d), intent(in) :: u
      type(amr_tree_state_3d), intent(inout) :: t
      logical, intent(out) :: good
      call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,u,t,comm,good)
    end subroutine
    subroutine sources(u,t,external,diffusive,face_flux,good)
      type(amr_tree_state_3d), intent(in) :: u,t
      type(amr_tree_state_3d), intent(inout) :: external,diffusive
      type(amr_tree_flux_3d), allocatable, intent(inout) :: face_flux(:)
      logical, intent(out) :: good
      integer :: peer
      external=u;diffusive=u;good=.true.
      do peer=1,size(tree%nodes)
        if(owners(peer)/=rank) cycle
        call compute_reactive_external_source_3d(u%fields(peer)%values,force,external%fields(peer)%values,good)
        if(.not.good) exit
        diffusive%fields(peer)%values=0
      end do
      call agree(good)
      if(.not.good) return
      if(do_transport) then
        stage_theta=1
        call compute_sdc_mpi_amr_tree_transport_3d(tree,owners,species,transport,u,t,spacing,dt,faces,transport_flags, &
          diffusive,face_flux,stage_theta,comm,good,inflow_state,inflow_temperature,wall_temperature,wall_velocity, &
          outlet_pressure,inflow_mach,inflow_pressure)
        if(.not.good) return
        theta=min(theta,stage_theta)
      end if
      good=.true.
    end subroutine
    subroutine agree(condition)
      logical, intent(inout) :: condition
      logical :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      condition=.false.
      if(ierr==MPI_SUCCESS) condition=result
    end subroutine
    subroutine inventory(fields,result,good)
      type(amr_tree_state_3d), intent(in) :: fields
      real(dp), intent(out) :: result(:)
      logical, intent(out) :: good
      real(dp) :: accum(size(result),2),term(size(result)),updated(size(result)),vol(size(tree%nodes)),ratio
      integer :: peer,parent,ii,jj,kk,dims(3),ierr
      good=.false.;accum=0;vol(1)=product(spacing)
      if(.not.ieee_is_finite(vol(1)).or.vol(1)<=0) return
      do peer=1,size(tree%nodes)
        if(peer>1) then
          parent=tree%nodes(peer)%parent;ratio=real(tree%nodes(peer)%patch%refinement_ratio,dp)
          vol(peer)=vol(parent)/ratio/ratio/ratio
          if(vol(peer)<=0) return
        end if
        if(owners(peer)==rank) then
          dims=tree%nodes(peer)%shape
          do kk=1,dims(3)
          do jj=1,dims(2)
          do ii=1,dims(1)
            if(covered(tree,peer,[ii,jj,kk])) cycle
            term=vol(peer)*fields%fields(peer)%values(:,ii,jj,kk)-accum(:,2)
            updated=accum(:,1)+term
            accum(:,2)=(updated-accum(:,1))-term;accum(:,1)=updated
          end do
          end do
          end do
        end if
        call MPI_Bcast(accum,size(accum),MPI_DOUBLE_PRECISION,owners(peer),comm,ierr)
        if(ierr/=MPI_SUCCESS) return
      end do
      if(.not.all(ieee_is_finite(accum))) return
      result=accum(:,1);good=.true.
    end subroutine
  end subroutine advance_mpi_amr_tree_sdc_3d

  ! R/2 - T/2 - H - T/2 - R/2, over one coarse interval, optionally with recursive flow/transport steps.
  ! Disabled operators are skipped. All stage outputs remain private until
  ! the complete composition succeeds, including the final chemistry solve.
  subroutine advance_mpi_amr_tree_split_3d(tree,owners,species,reactions,transport,state,temperature,spacing,dt,faces, &
      reconstruction,limiter,solver,chemistry,viscosity,conduction,diffusion,barodiffusion,rtol,atol, &
      boundary_change,maximum_correction,minimum_theta,comm,ok,chemistry_integrator, &
      inflow_state,inflow_temperature,wall_temperature,wall_velocity,subcycling,outlet_pressure,inflow_mach,inflow_pressure, &
      difmag)
    type(amr_patch_tree_3d), intent(in) :: tree
    integer, intent(in) :: owners(:)
    type(MPI_Comm), intent(in) :: comm
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
    integer :: rank,status,flags(6),low(6),high(6)
    ok=.false.
    valid=ieee_is_finite(dt).and.dt>0.and..not.(barodiffusion.and..not.diffusion)
    call agree(valid)
    if (.not.valid) return
    recursive_steps=.false.
    if (present(subcycling)) recursive_steps=subcycling
    flags=[merge(1,0,chemistry),merge(1,0,viscosity),merge(1,0,conduction),merge(1,0,diffusion), &
      merge(1,0,barodiffusion),merge(1,0,recursive_steps)]
    call MPI_Allreduce(flags,low,6,MPI_INTEGER,MPI_MIN,comm,status)
    if (status/=MPI_SUCCESS) return
    call MPI_Allreduce(flags,high,6,MPI_INTEGER,MPI_MAX,comm,status)
    if (status/=MPI_SUCCESS) return
    if (any(low/=high)) return
    if (reconstruction=='characteristic_plm'.or.reconstruction=='characteristic_ppm') transport_limiter=limiter
    call clone_mpi_amr_tree_inputs_3d(tree,owners,size(species),state,temperature,candidate,temps,comm,rank,valid)
    if (.not.valid) return
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
      call advance_mpi_amr_tree_subcycle_hydro_3d(tree,owners,species,candidate,temps,spacing,dt,faces, &
        reconstruction,limiter,solver,stage_budget,stage_max,comm,valid,inflow_state,inflow_temperature, &
        outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure,difmag=difmag)
    else
      call advance_mpi_amr_tree_hydro_3d(tree,owners,species,candidate,temps,spacing,dt,faces,reconstruction,limiter,solver, &
        stage_budget,stage_max,comm,valid,inflow_state,inflow_temperature, &
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
    call move_alloc(candidate%fields,state%fields)
    call move_alloc(temps%fields,temperature%fields)
    boundary_change=budget;maximum_correction=largest;minimum_theta=theta;ok=.true.
  contains
    subroutine reaction_half(success)
      logical, intent(out) :: success
      call advance_mpi_amr_tree_chemistry_3d(tree,owners,species,reactions,candidate,temps,.5_dp*dt,rtol,atol, &
        comm,success,chemistry_integrator)
    end subroutine
    subroutine transport_half(success)
      logical, intent(out) :: success
      stage_budget=0;stage_max=0;stage_theta=1
      if (recursive_steps) then
        call advance_mpi_amr_tree_subcycle_transport_3d(tree,owners,species,transport,candidate,temps,spacing, &
          .5_dp*dt,faces,viscosity,conduction,diffusion,barodiffusion,stage_budget,stage_max,stage_theta,comm,success, &
          inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter=transport_limiter, &
          outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
      else
        call advance_mpi_amr_tree_transport_3d(tree,owners,species,transport,candidate,temps,spacing,.5_dp*dt,faces, &
          viscosity,conduction,diffusion,barodiffusion,stage_budget,stage_max,stage_theta,comm,success, &
          inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter=transport_limiter, &
          outlet_pressure=outlet_pressure,inflow_mach=inflow_mach,inflow_pressure=inflow_pressure)
      end if
      if (.not.success) return
      budget=budget+stage_budget;largest=max(largest,stage_max);theta=min(theta,stage_theta)
    end subroutine
    subroutine agree(condition)
      logical, intent(inout) :: condition
      logical :: result
      integer :: ierr
      call MPI_Allreduce(condition,result,1,MPI_LOGICAL,MPI_LAND,comm,ierr)
      condition=.false.
      if (ierr==MPI_SUCCESS) condition=result
    end subroutine
  end subroutine advance_mpi_amr_tree_split_3d
end module mpi_amr_patch_tree_split_3d_mod
