module amr_patch_tree_hydro_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: reactive_nvar
  use reactive_3d_mod, only: recover_reactive_temperatures_3d,compute_reactive_cfl_timestep_3d, &
    compute_reactive_face_fluxes_3d,compute_reactive_plm_face_fluxes_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,average_down_amr_patch_tree_3d
  use amr_patch_tree_ghost_3d_mod, only: fill_amr_patch_tree_ghosts_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d,reflux_amr_patch_tree_3d
  use simulation_config_reactive_3d_mod, only: valid_reactive_boundaries_3d
  implicit none
  private
  public :: advance_amr_patch_tree_hydro_3d,compute_amr_patch_tree_cfl_3d
contains
  subroutine compute_amr_patch_tree_cfl_3d(tree,species,state,temperature,root_spacing,cfl,dt,ok)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(in) :: state,temperature
    real(dp), intent(in) :: root_spacing(3),cfl
    real(dp), intent(inout) :: dt
    logical, intent(out) :: ok
    real(dp), allocatable :: h(:,:)
    real(dp) :: candidate,local_dt
    integer :: node,n(3)
    logical :: valid
    ok=.false.
    call prepare(tree,species,state,temperature,root_spacing,h,valid)
    if (.not.valid) return
    candidate=huge(1.0_dp)
    do node=1,size(tree%nodes)
      n=tree%nodes(node)%shape
      call compute_reactive_cfl_timestep_3d(species,state%fields(node)%values,temperature%fields(node)%values(1,:,:,:), &
        n(1),n(2),n(3),h(1,node),h(2,node),h(3,node),cfl,local_dt,valid)
      if (.not.valid) return
      candidate=min(candidate,local_dt)
    end do
    dt=candidate;ok=.true.
  end subroutine
  subroutine advance_amr_patch_tree_hydro_3d(tree,species,state,temperature,root_spacing,dt,faces, &
      reconstruction,limiter,solver,boundary_change,maximum_correction,ok,inflow_state,inflow_temperature)
    type(amr_patch_tree_3d), intent(in) :: tree
    type(nasa7_species), intent(in) :: species(:)
    type(amr_tree_state_3d), intent(inout) :: state,temperature
    real(dp), intent(in) :: root_spacing(3),dt
    character(len=*), intent(in) :: faces(6),reconstruction,limiter,solver
    real(dp), intent(inout) :: boundary_change(:),maximum_correction
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:)
    type(amr_tree_state_3d) :: initial,candidate,temps,stage_temperature
    real(dp), allocatable :: h(:,:),inlet(:,:)
    real(dp) :: first_budget(size(boundary_change)),second_budget(size(boundary_change)),first_max,second_max
    integer :: node,nvar
    logical :: valid
    ok=.false.
    call prepare(tree,species,state,temperature,root_spacing,h,valid)
    if (.not.valid) return
    nvar=reactive_nvar(size(species))
    if (size(boundary_change)/=nvar.or..not.valid_reactive_boundaries_3d(faces)) return
    if (.not.ieee_is_finite(dt)) return
    if (dt<=0.or.(reconstruction/='pcm'.and.reconstruction/='characteristic_plm')) return
    if (limiter/='mc'.and.limiter/='minmod') return
    if (any(faces=='inflow')) then
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
    stage_temperature=temps
    call euler_stage(candidate,temps,second_budget,second_max,valid)
    if (.not.valid) return
    do node=1,size(tree%nodes)
      candidate%fields(node)%values=.5_dp*initial%fields(node)%values+.5_dp*candidate%fields(node)%values
    end do
    call average_down_amr_patch_tree_3d(tree,candidate,valid)
    if (.not.valid) return
    temps=stage_temperature
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
      real(dp), allocatable :: ghost(:,:,:,:),fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:)
      real(dp) :: correction(nvar),rhs(nvar)
      integer :: peer,n(3),dims(3),i,j,k
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
        if (reconstruction=='pcm') then
          call fill_amr_patch_tree_ghosts_3d(tree,packed,peer,faces,ghost,success,inflow_state=inlet)
        else
          call fill_amr_patch_tree_ghosts_3d(tree,packed,peer,faces,ghost,success,inflow_state=inlet, &
            coarse_limiter=limiter,species_count=size(species))
        end if
        if (.not.success) return
        allocate(fx(nvar,dims(1),dims(2),dims(3)),fy(nvar,dims(1),dims(2),dims(3)),fz(nvar,dims(1),dims(2),dims(3)))
        if (reconstruction=='pcm') then
          call compute_reactive_face_fluxes_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:), &
            dims(1),dims(2),dims(3),solver,fx,fy,fz,success)
        else
          call compute_reactive_plm_face_fluxes_3d(species,ghost(:nvar,:,:,:),ghost(nvar+1,:,:,:), &
            dims(1),dims(2),dims(3),limiter,solver,fx,fy,fz,success)
        end if
        if (.not.success) return
        do k=1,n(3)
          do j=1,n(2)
            do i=1,n(1)
              rhs=-(fx(:,i+2,j+2,k+2)-fx(:,i+1,j+2,k+2))/h(1,peer) &
                -(fy(:,i+2,j+2,k+2)-fy(:,i+2,j+1,k+2))/h(2,peer) &
                -(fz(:,i+2,j+2,k+2)-fz(:,i+2,j+2,k+1))/h(3,peer)
              updated%fields(peer)%values(:,i,j,k)=u%fields(peer)%values(:,i,j,k)+dt*rhs
            end do
          end do
        end do
        flux(peer)%x=dt*fx(:,2:n(1)+2,3:n(2)+2,3:n(3)+2)
        flux(peer)%y=dt*fy(:,3:n(1)+2,2:n(2)+2,3:n(3)+2)
        flux(peer)%z=dt*fz(:,3:n(1)+2,3:n(2)+2,2:n(3)+2)
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
  end subroutine advance_amr_patch_tree_hydro_3d
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
end module amr_patch_tree_hydro_3d_mod
