module les_amr_patch_tree_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use les_dynamic_smagorinsky_mod, only: compute_dynamic_les_rhs_ghosted_3d
  use les_filter_3d_mod, only: les_filter_weights
  use les_smagorinsky_3d_mod, only: compute_les_ghosted_rhs_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,average_down_amr_patch_tree_3d, &
    amr_tree_integer_layout_3d
  use amr_patch_tree_ghost_3d_mod, only: fill_amr_patch_tree_ghosts_3d,apply_amr_tree_pressure_ghosts_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d,reflux_amr_patch_tree_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  implicit none
  private
  public :: advance_amr_patch_tree_les_3d,compute_amr_patch_tree_les_timestep_3d
  public :: advance_amr_patch_tree_les_level_3d
contains
  subroutine compute_amr_patch_tree_les_timestep_3d(tree,species,state,temperature,root_spacing,cfl,faces, &
      cs,prt,sct,dt,maximum_diffusivity,ok,inflow_state,inflow_temperature,wall_temperature,wall_velocity, &
      coarse_limiter,outlet_pressure,model,ci,test_model,test_ratio,grid_ratio,harmonic_mean,subcycle,inflow_mach,inflow_pressure)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: root_spacing(3),cfl,cs,prt,sct
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(inout) :: dt,maximum_diffusivity
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_pressure(6)
    real(dp) :: mach(6)
    real(dp), optional, intent(in) :: outlet_pressure(6),ci
    character(len=*), optional, intent(in) :: coarse_limiter,model
    character(len=*), optional, intent(in) :: test_model
    integer, optional, intent(in) :: test_ratio,grid_ratio
    logical, optional, intent(in) :: harmonic_mean,subcycle
    character(len=32) :: filter
    integer :: filter_ratio,width
    logical :: dynamic
    type(amr_tree_state_3d) :: averaged,packed
    real(dp), allocatable :: h(:,:),inlet(:,:),ghost(:,:,:,:),rhs(:,:,:,:),recovered(:,:,:)
    integer(int64), allocatable :: origin(:,:),domain(:,:),counts(:)
    integer, allocatable :: ratios(:)
    real(dp) :: candidate,diffusivity,largest,rate
    character(len=32) :: local_faces(6)
    integer :: peer,n(3),nvar,axis,level
    integer(int64) :: factor
    logical :: valid
    ok=.false.
    call dynamic_controls(model,faces,test_model,test_ratio,grid_ratio,harmonic_mean, &
      dynamic,filter,filter_ratio,width,valid)
    if(.not.valid) return
    call prepare(tree,species,state,temperature,root_spacing,h,valid)
    if (.not.valid) return
    if (.not.ieee_is_finite(cfl)) return
    if (cfl<=0.or.cfl>.2_dp.or..not.valid_reactive_boundaries_3d(faces)) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    if (.not.valid) return
    allocate(counts(size(ratios)));counts=1
    if(present(subcycle)) then
      if(subcycle) then
        do level=2,size(ratios)
          factor=int(ratios(level),int64)**2
          if(counts(level-1)>huge(1_int64)/factor) return
          counts(level)=counts(level-1)*factor
        end do
      end if
    end if
    nvar=reactive_nvar(size(species))
    mach=-1
    if(present(inflow_mach)) mach=inflow_mach
    if(present(inflow_temperature)) then
      if(size(inflow_temperature)/=6) return
    end if
    if (any(faces=='inflow'.and.mach<0)) then
      if (.not.present(inflow_state).or..not.present(inflow_temperature)) return
      if (any(shape(inflow_state)/=[nvar,6]).or.size(inflow_temperature)/=6) return
      if (.not.all(ieee_is_finite(inflow_state)).or..not.all(ieee_is_finite(inflow_temperature))) return
      allocate(inlet(nvar+1,6));inlet(:nvar,:)=inflow_state;inlet(nvar+1,:)=inflow_temperature
    end if
    averaged=state
    call average_down_amr_patch_tree_3d(tree,averaged,valid)
    if (.not.valid) return
    allocate(packed%fields(size(tree%nodes)))
    do peer=1,size(tree%nodes)
      n=tree%nodes(peer)%shape
      allocate(recovered(n(1),n(2),n(3)))
      call recover_reactive_temperatures_3d(species,averaged%fields(peer)%values, &
        temperature%fields(peer)%values(1,:,:,:),n(1),n(2),n(3),recovered,valid)
      if (.not.valid) return
      allocate(packed%fields(peer)%values(nvar+1,n(1),n(2),n(3)))
      packed%fields(peer)%values(:nvar,:,:,:)=averaged%fields(peer)%values
      packed%fields(peer)%values(nvar+1,:,:,:)=recovered
      deallocate(recovered)
    end do
    candidate=1e6_dp;largest=0
    do peer=1,size(tree%nodes)
      n=tree%nodes(peer)%shape
      call fill_amr_patch_tree_ghosts_3d(tree,packed,peer,faces,ghost,valid,viscous_wall=.true., &
        inflow_state=inlet,wall_velocity=wall_velocity,coarse_limiter=coarse_limiter, &
        species_count=size(species),ghost_width=width, &
        species=species,inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure)
      if (.not.valid) return
      call apply_amr_tree_pressure_ghosts_3d(tree,species,peer,faces,outlet_pressure,ghost,valid,ghost_width=width)
      if (.not.valid) return
      local_faces='outflow';level=tree%nodes(peer)%level+1
      do axis=1,3
        if (origin(axis,peer)==0) local_faces(2*axis-1)=faces(2*axis-1)
        if (origin(axis,peer)+int(n(axis),int64)==domain(axis,level)) local_faces(2*axis)=faces(2*axis)
      end do
      allocate(rhs(nvar,n(1),n(2),n(3)))
      if(dynamic) then
        block
          real(dp), allocatable :: fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:)
          allocate(fx(nvar,n(1)+1,n(2),n(3)),fy(nvar,n(1),n(2)+1,n(3)),fz(nvar,n(1),n(2),n(3)+1))
          call compute_dynamic_les_rhs_ghosted_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:), &
            h(:,peer),filter,filter_ratio,width,rhs,fx,fy,fz,valid,grid_ratio,harmonic_mean,diffusivity, &
              boundary_faces=local_faces,wall_temperature=wall_temperature)
        end block
      else
      call compute_les_ghosted_rhs_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:),n,h(:,peer), &
        cs,prt,sct,rhs,diffusivity,valid,boundary_faces=local_faces,wall_temperature=wall_temperature,model=model,ci=ci)
      end if
      if (.not.valid) return
      rate=diffusivity*sum(1._dp/h(:,peer)**2)/real(counts(level),dp)
      if (rate>cfl/candidate) candidate=cfl/rate
      largest=max(largest,diffusivity)
      deallocate(rhs)
    end do
    dt=candidate;maximum_diffusivity=largest;ok=.true.
  end subroutine

  subroutine advance_amr_patch_tree_les_3d(tree,species,state,temperature,root_spacing,dt,faces, &
      cs,prt,sct,boundary_change,maximum_correction,ok, &
      inflow_state,inflow_temperature,wall_temperature,wall_velocity,coarse_limiter,outlet_pressure,model,ci, &
      test_model,test_ratio,grid_ratio,harmonic_mean,inflow_mach,inflow_pressure)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: root_spacing(3),dt
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: cs,prt,sct
    character(len=*), optional, intent(in) :: model
    real(dp), optional, intent(in) :: ci
    real(dp), intent(inout) :: boundary_change(:),maximum_correction
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_pressure(6)
    real(dp) :: mach(6)
    real(dp), optional, intent(in) :: outlet_pressure(6)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    character(len=*), optional, intent(in) :: coarse_limiter
    character(len=*), optional, intent(in) :: test_model
    integer, optional, intent(in) :: test_ratio,grid_ratio
    logical, optional, intent(in) :: harmonic_mean
    character(len=32) :: filter
    integer :: filter_ratio,width
    logical :: dynamic
    type(amr_tree_state_3d) :: initial,candidate,temps,first_temperature
    real(dp), allocatable :: h(:,:),inlet(:,:)
    real(dp) :: first_budget(size(boundary_change)),second_budget(size(boundary_change)), &
      first_max,second_max
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    integer :: node,nvar
    logical :: valid
    ok=.false.
    call dynamic_controls(model,faces,test_model,test_ratio,grid_ratio,harmonic_mean, &
      dynamic,filter,filter_ratio,width,valid)
    if(.not.valid) return
    call prepare(tree,species,state,temperature,root_spacing,h,valid)
    if (.not.valid) return
    nvar=reactive_nvar(size(species))
    if (size(boundary_change)/=nvar.or..not.valid_reactive_boundaries_3d(faces)) return
    if (.not.ieee_is_finite(dt)) return
    if (dt<=0.or.dt>1e6_dp) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    if (.not.valid) return
    mach=-1
    if(present(inflow_mach)) mach=inflow_mach
    if(present(inflow_temperature)) then
      if(size(inflow_temperature)/=6) return
    end if
    if (any(faces=='inflow'.and.mach<0)) then
      if (.not.present(inflow_state).or..not.present(inflow_temperature)) return
      if (any(shape(inflow_state)/=[nvar,6]).or.size(inflow_temperature)/=6) return
      if (.not.all(ieee_is_finite(inflow_state)).or..not.all(ieee_is_finite(inflow_temperature))) return
      allocate(inlet(nvar+1,6));inlet(:nvar,:)=inflow_state;inlet(nvar+1,:)=inflow_temperature
    end if
    initial=state;temps=temperature
    call average_down_amr_patch_tree_3d(tree,initial,valid)
    if (.not.valid) return
    call recover_tree(initial,temps,valid)
    if (.not.valid) return
    candidate=initial
    call euler_stage(candidate,temps,first_budget,first_max,valid)
    if (.not.valid) return
    first_temperature=temps
    call euler_stage(candidate,temps,second_budget,second_max,valid)
    if (.not.valid) return
    do node=1,size(tree%nodes)
      candidate%fields(node)%values=.5_dp*initial%fields(node)%values+.5_dp*candidate%fields(node)%values
    end do
    call average_down_amr_patch_tree_3d(tree,candidate,valid)
    if (.not.valid) return
    temps=first_temperature
    call recover_tree(candidate,temps,valid)
    if (.not.valid) return
    state=candidate;temperature=temps
    boundary_change=.5_dp*(first_budget+second_budget)
    maximum_correction=max(first_max,second_max);ok=.true.
  contains
    subroutine recover_tree(u,t,success)
      type(amr_tree_state_3d), intent(in) :: u
      type(amr_tree_state_3d), intent(inout) :: t
      logical, intent(out) :: success
      real(dp), allocatable :: recovered(:,:,:)
      integer :: peer,n(3)
      do peer=1,size(tree%nodes)
        n=tree%nodes(peer)%shape
        allocate(recovered(n(1),n(2),n(3)))
        call recover_reactive_temperatures_3d(species,u%fields(peer)%values,t%fields(peer)%values(1,:,:,:), &
          n(1),n(2),n(3),recovered,success)
        if (.not.success) return
        t%fields(peer)%values(1,:,:,:)=recovered
        deallocate(recovered)
      end do
      success=.true.
    end subroutine
    subroutine euler_stage(u,t,budget,largest,success)
      type(amr_tree_state_3d), intent(inout) :: u,t
      real(dp), intent(out) :: budget(:),largest
      logical, intent(out) :: success
      type(amr_tree_state_3d) :: packed,updated
      type(amr_tree_flux_3d), allocatable :: flux(:)
      real(dp), allocatable :: ghost(:,:,:,:),fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:),les_rhs(:,:,:,:)
      real(dp) :: correction(nvar),rhs(nvar),diffusivity
      character(len=32) :: local_faces(6)
      integer :: peer,n(3),dims(3),i,j,k,axis,level
      success=.false.;budget=0;largest=0;correction=0
      allocate(packed%fields(size(tree%nodes)),flux(size(tree%nodes)))
      do peer=1,size(tree%nodes)
        n=tree%nodes(peer)%shape
        allocate(packed%fields(peer)%values(nvar+1,n(1),n(2),n(3)))
        packed%fields(peer)%values(:nvar,:,:,:)=u%fields(peer)%values
        packed%fields(peer)%values(nvar+1,:,:,:)=t%fields(peer)%values(1,:,:,:)
      end do
      updated=u
      do peer=1,size(tree%nodes)
        n=tree%nodes(peer)%shape;dims=n+4
        call fill_amr_patch_tree_ghosts_3d(tree,packed,peer,faces,ghost,success, &
          viscous_wall=.true.,inflow_state=inlet,wall_velocity=wall_velocity, &
          coarse_limiter=coarse_limiter,species_count=size(species),ghost_width=width, &
        species=species,inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure)
        if (.not.success) return
        call apply_amr_tree_pressure_ghosts_3d(tree,species,peer,faces,outlet_pressure,ghost,success,ghost_width=width)
        if (.not.success) return
        allocate(fx(nvar,n(1)+1,n(2),n(3)),fy(nvar,n(1),n(2)+1,n(3)),fz(nvar,n(1),n(2),n(3)+1))
        local_faces='outflow';level=tree%nodes(peer)%level+1
        do axis=1,3
          if (origin(axis,peer)==0) local_faces(2*axis-1)=faces(2*axis-1)
          if (origin(axis,peer)+int(n(axis),int64)==domain(axis,level)) local_faces(2*axis)=faces(2*axis)
        end do
        allocate(les_rhs(nvar,n(1),n(2),n(3)))
        if(dynamic) then
          call compute_dynamic_les_rhs_ghosted_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:), &
            h(:,peer),filter,filter_ratio,width,les_rhs,fx,fy,fz,success,grid_ratio,harmonic_mean,diffusivity, &
              boundary_faces=local_faces,wall_temperature=wall_temperature)
        else
        call compute_les_ghosted_rhs_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:),n,h(:,peer), &
          cs,prt,sct,les_rhs,diffusivity,success,boundary_faces=local_faces, &
          wall_temperature=wall_temperature,model=model,ci=ci,fx=fx,fy=fy,fz=fz)
        end if
        if (.not.success) return
        success=dt*diffusivity*sum(1._dp/h(:,peer)**2)<=.2_dp
        if (.not.success) return
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
        flux(peer)%x=dt*fx
        flux(peer)%y=dt*fy
        flux(peer)%z=dt*fz
        deallocate(fx,fy,fz)
      end do
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
      call reflux_amr_patch_tree_3d(tree,updated,flux,root_spacing,faces,correction,largest,success)
      if (.not.success) return
      budget=budget+correction
      call recover_tree(updated,t,success)
      if (.not.success) return
      u=updated
    end subroutine
  end subroutine advance_amr_patch_tree_les_3d
  ! Advances only depth, returning SSPRK2 time-integrated faces for later reflux.
  ! Coarse endpoint fields pack conserved state followed by a temperature guess.
  ! The recursive caller owns synchronization/restriction after child intervals.
  subroutine advance_amr_patch_tree_les_level_3d(tree,species,state,temperature,depth,root_spacing,dt,faces, &
      cs,prt,sct,integrated_flux,maximum_diffusivity,ok, &
      coarse_start,coarse_end,fractions_start,fractions_end,inflow_state,inflow_temperature, &
      wall_temperature,wall_velocity,coarse_limiter,outlet_pressure,model,ci, &
      test_model,test_ratio,grid_ratio,harmonic_mean,inflow_mach,inflow_pressure)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    integer, intent(in) :: depth
    real(dp), intent(in) :: root_spacing(3),dt
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: cs,prt,sct
    character(len=*), optional, intent(in) :: model,test_model
    real(dp), optional, intent(in) :: ci
    integer, optional, intent(in) :: test_ratio,grid_ratio
    logical, optional, intent(in) :: harmonic_mean
    character(len=32) :: filter
    integer :: filter_ratio,width
    logical :: dynamic
    real(dp), intent(inout) :: maximum_diffusivity
    real(dp), optional, intent(in) :: wall_temperature(:),wall_velocity(:,:)
    character(len=*), optional, intent(in) :: coarse_limiter
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer, allocatable :: ratios(:)
    character(len=32) :: local_faces(6)
    real(dp) :: largest,diffusivity
    type(amr_tree_flux_3d), allocatable, intent(inout) :: integrated_flux(:)
    real(dp), optional, intent(in) :: inflow_mach(6),inflow_pressure(6)
    real(dp) :: mach(6)
    real(dp), optional, intent(in) :: outlet_pressure(6)
    logical, intent(out) :: ok
    ! Endpoints pack conserved components followed by a temperature guess.
    type(amr_tree_state_3d), optional, intent(in) :: coarse_start,coarse_end
    real(dp), optional, intent(in) :: fractions_start(:),fractions_end(:),inflow_state(:,:),inflow_temperature(:)
    type(amr_tree_state_3d) :: initial,candidate,temps,packed,updated,first_t
    type(amr_tree_flux_3d), allocatable :: flux(:)
    real(dp), allocatable :: h(:,:),inlet(:,:),ghost(:,:,:,:),les_rhs(:,:,:,:),fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:),recovered(:,:,:)
    real(dp) :: rhs(reactive_nvar(size(species)))
    integer :: peer,n(3),dims(3),nvar,stage,i,j,k,endpoint,axis,level
    logical :: valid
    ok=.false.
    call dynamic_controls(model,faces,test_model,test_ratio,grid_ratio,harmonic_mean, &
      dynamic,filter,filter_ratio,width,valid)
    if(.not.valid) return
    call prepare(tree,species,state,temperature,root_spacing,h,valid)
    if (.not.valid) return
    if (depth<0.or.depth>=tree%level_count()) return
    if (.not.valid_reactive_boundaries_3d(faces).or..not.ieee_is_finite(dt)) return
    if (dt<=0.or.dt>1e6_dp) return
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,valid)
    if (.not.valid) return
    largest=0
    nvar=reactive_nvar(size(species))
    if (depth>0) then
      if (.not.present(coarse_start).or..not.present(coarse_end).or. &
          .not.present(fractions_start).or..not.present(fractions_end)) return
    end if
    if (present(coarse_start).neqv.present(coarse_end)) return
    if (present(coarse_end).neqv.present(fractions_start)) return
    if (present(fractions_start).neqv.present(fractions_end)) return
    if (present(coarse_start)) then
      if (size(fractions_start)/=depth.or.size(fractions_end)/=depth) return
      if (.not.all(ieee_is_finite(fractions_start)).or..not.all(ieee_is_finite(fractions_end))) return
      if (any(fractions_start<0).or.any(fractions_end>1).or.any(fractions_end<fractions_start)) return
      if (.not.allocated(coarse_start%fields).or..not.allocated(coarse_end%fields)) return
      if (size(coarse_start%fields)/=size(tree%nodes).or.size(coarse_end%fields)/=size(tree%nodes)) return
      do peer=1,size(tree%nodes)
        if (tree%nodes(peer)%level>=depth) cycle
        n=tree%nodes(peer)%shape
        if (.not.allocated(coarse_start%fields(peer)%values).or..not.allocated(coarse_end%fields(peer)%values)) return
        if (any(shape(coarse_start%fields(peer)%values)/=[nvar+1,n]).or. &
            any(shape(coarse_end%fields(peer)%values)/=[nvar+1,n])) return
        if (.not.all(ieee_is_finite(coarse_start%fields(peer)%values)).or. &
            .not.all(ieee_is_finite(coarse_end%fields(peer)%values))) return
      end do
    end if
    mach=-1
    if(present(inflow_mach)) mach=inflow_mach
    if(present(inflow_temperature)) then
      if(size(inflow_temperature)/=6) return
    end if
    if (any(faces=='inflow'.and.mach<0)) then
      if (.not.present(inflow_state).or..not.present(inflow_temperature)) return
      if (any(shape(inflow_state)/=[nvar,6]).or.size(inflow_temperature)/=6) return
      if (.not.all(ieee_is_finite(inflow_state)).or..not.all(ieee_is_finite(inflow_temperature))) return
      allocate(inlet(nvar+1,6));inlet(:nvar,:)=inflow_state;inlet(nvar+1,:)=inflow_temperature
    end if
    initial=state;candidate=state;temps=temperature
    allocate(flux(size(tree%nodes)),packed%fields(size(tree%nodes)))
    do stage=1,2
      do peer=1,size(tree%nodes)
        n=tree%nodes(peer)%shape
        if (.not.allocated(packed%fields(peer)%values)) allocate(packed%fields(peer)%values(nvar+1,n(1),n(2),n(3)))
        packed%fields(peer)%values(:nvar,:,:,:)=candidate%fields(peer)%values
        packed%fields(peer)%values(nvar+1,:,:,:)=temps%fields(peer)%values(1,:,:,:)
        if (tree%nodes(peer)%level<depth) packed%fields(peer)%values=coarse_start%fields(peer)%values
      end do
      updated=candidate
      do peer=1,size(tree%nodes)
        if (tree%nodes(peer)%level/=depth) cycle
        n=tree%nodes(peer)%shape;dims=n+4
        endpoint=stage
        call stage_ghost(peer,valid)
        if (.not.valid) return
        allocate(fx(nvar,n(1)+1,n(2),n(3)),fy(nvar,n(1),n(2)+1,n(3)),fz(nvar,n(1),n(2),n(3)+1))
        local_faces='outflow';level=tree%nodes(peer)%level+1
        do axis=1,3
          if (origin(axis,peer)==0) local_faces(2*axis-1)=faces(2*axis-1)
          if (origin(axis,peer)+int(n(axis),int64)==domain(axis,level)) local_faces(2*axis)=faces(2*axis)
        end do
        allocate(les_rhs(nvar,n(1),n(2),n(3)))
        if(dynamic) then
          call compute_dynamic_les_rhs_ghosted_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:), &
            h(:,peer),filter,filter_ratio,width,les_rhs,fx,fy,fz,valid,grid_ratio,harmonic_mean,diffusivity, &
            boundary_faces=local_faces,wall_temperature=wall_temperature)
        else
          call compute_les_ghosted_rhs_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:),n,h(:,peer), &
            cs,prt,sct,les_rhs,diffusivity,valid,boundary_faces=local_faces, &
            wall_temperature=wall_temperature,model=model,ci=ci,fx=fx,fy=fy,fz=fz)
        end if
        if(.not.valid) return
        if(dt*diffusivity*sum(1._dp/h(:,peer)**2)>.2_dp) return
        largest=max(largest,diffusivity)
        deallocate(les_rhs)
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              rhs=-(fx(:,i+1,j,k)-fx(:,i,j,k))/h(1,peer) &
                -(fy(:,i,j+1,k)-fy(:,i,j,k))/h(2,peer) &
                -(fz(:,i,j,k+1)-fz(:,i,j,k))/h(3,peer)
              updated%fields(peer)%values(:,i,j,k)=candidate%fields(peer)%values(:,i,j,k)+dt*rhs
            end do
          end do
        end do
        if (stage==1) then
          flux(peer)%x=.5_dp*dt*fx
          flux(peer)%y=.5_dp*dt*fy
          flux(peer)%z=.5_dp*dt*fz
        else
          flux(peer)%x=flux(peer)%x+.5_dp*dt*fx
          flux(peer)%y=flux(peer)%y+.5_dp*dt*fy
          flux(peer)%z=flux(peer)%z+.5_dp*dt*fz
        end if
        deallocate(fx,fy,fz)
      end do
      do peer=1,size(tree%nodes)
        if (tree%nodes(peer)%level/=depth) cycle
        n=tree%nodes(peer)%shape;allocate(recovered(n(1),n(2),n(3)))
        call recover_reactive_temperatures_3d(species,updated%fields(peer)%values,temps%fields(peer)%values(1,:,:,:), &
          n(1),n(2),n(3),recovered,valid)
        if (.not.valid) return
        temps%fields(peer)%values(1,:,:,:)=recovered;deallocate(recovered)
      end do
      candidate=updated
      if(stage==1) first_t=temps

    end do
    do peer=1,size(tree%nodes)
      if (tree%nodes(peer)%level/=depth) cycle
      candidate%fields(peer)%values=.5_dp*initial%fields(peer)%values+.5_dp*candidate%fields(peer)%values
      temps%fields(peer)%values=first_t%fields(peer)%values
      n=tree%nodes(peer)%shape;allocate(recovered(n(1),n(2),n(3)))
      call recover_reactive_temperatures_3d(species,candidate%fields(peer)%values,temps%fields(peer)%values(1,:,:,:), &
        n(1),n(2),n(3),recovered,valid)
      if (.not.valid) return
      temps%fields(peer)%values(1,:,:,:)=recovered;deallocate(recovered)
    end do
    do peer=1,size(tree%nodes)
      if (tree%nodes(peer)%level/=depth) cycle
      if (.not.all(ieee_is_finite(flux(peer)%x)).or..not.all(ieee_is_finite(flux(peer)%y)).or. &
          .not.all(ieee_is_finite(flux(peer)%z))) return
    end do
    state=candidate;temperature=temps;maximum_diffusivity=largest;call move_alloc(flux,integrated_flux);ok=.true.
  contains
    subroutine stage_ghost(peer,success)
      integer, intent(in) :: peer
      logical, intent(out) :: success
      real(dp), allocatable :: fractions(:)
      if (present(fractions_start)) then
        if (endpoint==1) then
          fractions=fractions_start
        else
          fractions=fractions_end
        end if
      end if
      call fill_amr_patch_tree_ghosts_3d(tree,packed,peer,faces,ghost,success,inflow_state=inlet, &
        coarse_future=coarse_end,coarse_fractions=fractions,coarse_limiter=coarse_limiter,species_count=size(species), &
        viscous_wall=.true.,wall_velocity=wall_velocity,ghost_width=width, &
        species=species,inflow_mach=inflow_mach,inflow_temperature=inflow_temperature,inflow_pressure=inflow_pressure)
      if (.not.success) return
      call apply_amr_tree_pressure_ghosts_3d(tree,species,peer,faces,outlet_pressure,ghost,success,ghost_width=width)
    end subroutine
  end subroutine advance_amr_patch_tree_les_level_3d

  subroutine dynamic_controls(model,faces,test_model,test_ratio,grid_ratio,harmonic_mean, &
      dynamic,filter,ratio,width,ok)
    character(len=*), optional, intent(in) :: model,test_model
    character(len=*), intent(in) :: faces(6)
    integer, optional, intent(in) :: test_ratio,grid_ratio
    logical, optional, intent(in) :: harmonic_mean
    logical, intent(out) :: dynamic,ok
    character(len=*), intent(out) :: filter
    integer, intent(out) :: ratio,width
    real(dp), allocatable :: weights(:)
    dynamic=.false.;ok=.false.;filter='box';ratio=2;width=2
    if(present(model)) dynamic=model=='dynamic'
    if(.not.dynamic) then
      if(present(test_model).or.present(test_ratio).or.present(grid_ratio).or.present(harmonic_mean)) return
      ok=.true.;return
    end if
    if(any(faces/='periodic'.and.faces/='inflow'.and.faces/='outflow'.and. &
      faces/='slip_wall'.and.faces/='no_slip_wall'.and.faces/='isothermal_wall')) return
    if(present(test_model)) filter=test_model
    if(present(test_ratio)) ratio=test_ratio
    if(present(grid_ratio)) then
      if(grid_ratio<1.or.grid_ratio>64) return
    end if
    call les_filter_weights(filter,ratio,weights,ok)
    if(.not.ok) return
    width=(size(weights)-1)/2+5
  end subroutine
  subroutine prepare(tree,species,state,temperature,root_spacing,h,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: root_spacing(3)
    real(dp), allocatable, intent(out) :: h(:,:)
    logical, intent(out) :: ok
    integer :: node,parent
    ok=.false.
    if (.not.state%is_valid(tree).or..not.temperature%is_valid(tree)) return
    if (size(state%fields(1)%values,1)/=reactive_nvar(size(species))) return
    if (size(temperature%fields(1)%values,1)/=1) return
    if (.not.all(ieee_is_finite(root_spacing))) return
    if (any(root_spacing<=0)) return
    allocate(h(3,size(tree%nodes)));h(:,1)=root_spacing
    do node=1,size(tree%nodes)
      if (any(temperature%fields(node)%values<=0)) return
      if (node==1) cycle
      parent=tree%nodes(node)%parent
      h(:,node)=h(:,parent)/real(tree%nodes(node)%patch%refinement_ratio,dp)
      if (any(h(:,node)<=0)) return
    end do
    ok=.true.
  end subroutine
end module les_amr_patch_tree_mod
