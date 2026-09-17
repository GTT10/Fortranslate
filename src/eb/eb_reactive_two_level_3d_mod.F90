module eb_reactive_two_level_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use gas_transport_mod, only: gas_transport_species
  use eb_reactive_transport_3d_mod, only: advance_reactive_eb_transport_patch_euler_3d
  use elementary_kinetics_mod, only: elementary_reaction
  use reactive_3d_mod, only: advance_reactive_chemistry_3d
  use reactive_1d_mod, only: reactive_nvar
  use eb_geometry_3d_mod, only: eb_geometry_3d,coarsen_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_transfer_3d_mod, only: restrict_reactive_eb_state_3d,fill_reactive_eb_patch_halo_3d
  use eb_reactive_hydro_3d_mod, only: advance_reactive_eb_euler_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d,restrict_eb_integrated_flux_3d,reflux_reactive_eb_state_3d
  use eb_reactive_wall_flux_3d_mod, only: reactive_eb_slip_wall_source_3d
  implicit none
  private
  public :: advance_reactive_eb_two_level_euler_3d,advance_reactive_eb_two_level_strang_3d
  public :: advance_reactive_eb_two_level_transport_3d,advance_reactive_eb_two_level_full_3d
contains
  ! Symmetric reaction/hydro composition; the Euler hydro stage remains first
  ! order. Only exposed coarse cells and active fine cells own chemical solves.
  subroutine advance_reactive_eb_two_level_strang_3d(species,reactions,coarse,fine,lower,ratio,solver,dt, &
      chemistry_enabled,rtol,atol,coarse_state,coarse_temperature,fine_state,fine_temperature, &
      boundary_change,wall_change,chemistry_change,ok,chemistry_integrator)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(eb_geometry_3d), intent(in) :: coarse,fine
    integer, intent(in) :: lower(3),ratio
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt,rtol,atol
    logical, intent(in) :: chemistry_enabled
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:)
    real(dp), intent(inout) :: fine_state(:,:,:,:),fine_temperature(:,:,:)
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemistry_change(:)
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: chemistry_integrator
    call advance_two_level_split(species,reactions,coarse,fine,lower,ratio,solver,dt,chemistry_enabled,rtol,atol, &
      coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall_change,chemistry_change, &
      ok,chemistry_integrator)
  end subroutine advance_reactive_eb_two_level_strang_3d

  subroutine advance_reactive_eb_two_level_full_3d(species,reactions,transport,coarse,fine,lower,ratio,solver,dt, &
      chemistry_enabled,rtol,atol,transport_enabled,viscosity,conduction,diffusion,barodiffusion, &
      coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall_change,chemistry_change, &
      minimum_theta,ok,chemistry_integrator)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_geometry_3d), intent(in) :: coarse,fine
    integer, intent(in) :: lower(3),ratio
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt,rtol,atol
    logical, intent(in) :: chemistry_enabled,transport_enabled,viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:),fine_state(:,:,:,:),fine_temperature(:,:,:)
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemistry_change(:),minimum_theta
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: chemistry_integrator
    if (transport_enabled) then
      call advance_two_level_split(species,reactions,coarse,fine,lower,ratio,solver,dt,chemistry_enabled,rtol,atol, &
        coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall_change,chemistry_change, &
        ok,chemistry_integrator,transport,viscosity,conduction,diffusion,barodiffusion,minimum_theta)
    else
      call advance_two_level_split(species,reactions,coarse,fine,lower,ratio,solver,dt,chemistry_enabled,rtol,atol, &
        coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall_change,chemistry_change, &
        ok,chemistry_integrator)
      if (ok) minimum_theta=1._dp
    end if
  end subroutine advance_reactive_eb_two_level_full_3d

  subroutine advance_two_level_split(species,reactions,coarse,fine,lower,ratio,solver,dt, &
      chemistry_enabled,rtol,atol,coarse_state,coarse_temperature,fine_state,fine_temperature, &
      boundary_change,wall_change,chemistry_change,ok,chemistry_integrator, &
      transport,viscosity,conduction,diffusion,barodiffusion,minimum_theta)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(eb_geometry_3d), intent(in) :: coarse,fine
    integer, intent(in) :: lower(3),ratio
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt,rtol,atol
    logical, intent(in) :: chemistry_enabled
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:)
    real(dp), intent(inout) :: fine_state(:,:,:,:),fine_temperature(:,:,:)
    real(dp), intent(inout) :: boundary_change(:),wall_change(:),chemistry_change(:)
    logical, intent(out) :: ok
    character(len=*), optional, intent(in) :: chemistry_integrator
    type(gas_transport_species), optional, intent(in) :: transport(:)
    logical, optional, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), optional, intent(inout) :: minimum_theta
    real(dp), allocatable :: transport_boundary(:),transport_sum(:)
    real(dp) :: stage_theta,theta
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:),child(:,:,:,:),child_t(:,:,:)
    real(dp), allocatable :: boundary(:),wall(:),chemical(:),before(:)
    logical, allocatable :: active_coarse(:,:,:),active_fine(:,:,:)
    integer :: nvar,n(3),nf(3),upper(3)
    logical :: valid
    ok=.false.;nvar=reactive_nvar(size(species))
    if (size(chemistry_change)/=nvar) return
    if (.not.chemistry_enabled.and..not.present(transport)) then
      call advance_reactive_eb_two_level_euler_3d(species,coarse,fine,lower,ratio,solver,dt, &
        coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall_change,ok)
      if (ok) chemistry_change=0._dp
      return
    end if
    if (.not.ieee_is_finite(dt)) return
    if (dt<=0._dp.or.ratio<1) return
    if (chemistry_enabled) then
      if (.not.all(ieee_is_finite([rtol,atol]))) return
      if (min(rtol,atol)<=0._dp.or.size(reactions)<1) return
    end if
    if (.not.coarse%is_valid().or..not.fine%is_valid()) return
    n=[coarse%nx,coarse%ny,coarse%nz];nf=[fine%nx,fine%ny,fine%nz]
    if (any(mod(nf,ratio)/=0).or.any(lower<1).or.any(lower>n)) return
    if (any(nf/ratio>n-lower+1)) return
    upper=lower+nf/ratio-1
    if (any(shape(coarse_state)/=[nvar,n]).or.any(shape(coarse_temperature)/=n).or. &
      any(shape(fine_state)/=[nvar,nf]).or.any(shape(fine_temperature)/=nf)) return
    if (size(boundary_change)/=nvar.or.size(wall_change)/=nvar) return
    allocate(candidate,source=coarse_state);allocate(temps,source=coarse_temperature)
    allocate(child,source=fine_state);allocate(child_t,source=fine_temperature)
    allocate(active_coarse(n(1),n(2),n(3)),active_fine(nf(1),nf(2),nf(3)))
    active_coarse=coarse%volume_fraction>0._dp;active_fine=fine%volume_fraction>0._dp
    active_coarse(lower(1):upper(1),lower(2):upper(2),lower(3):upper(3))=.false.
    allocate(boundary(nvar),wall(nvar),chemical(nvar),before(nvar));chemical=0._dp
    if (chemistry_enabled) then
      call reaction_half(valid)
      if (.not.valid) return
    end if
    allocate(transport_boundary(nvar),transport_sum(nvar));transport_sum=0._dp;theta=1._dp
    if (present(transport)) then
      call transport_half(valid)
      if (.not.valid) return
    end if
    call advance_reactive_eb_two_level_euler_3d(species,coarse,fine,lower,ratio,solver,dt, &
      candidate,temps,child,child_t,boundary,wall,valid)
    if (.not.valid) return
    if (present(transport)) then
      call transport_half(valid)
      if (.not.valid) return
    end if
    if (chemistry_enabled) then
      call reaction_half(valid)
      if (.not.valid) return
    end if
    if (present(transport)) boundary=boundary+transport_sum
    if (.not.all(ieee_is_finite(chemical)).or..not.all(ieee_is_finite(boundary))) return
    coarse_state=candidate;coarse_temperature=temps;fine_state=child;fine_temperature=child_t
    boundary_change=boundary;wall_change=wall;chemistry_change=chemical
    if (present(minimum_theta)) minimum_theta=theta
    ok=.true.
  contains
    subroutine transport_half(success)
      logical, intent(out) :: success
      call advance_reactive_eb_two_level_transport_3d(species,transport,coarse,fine,lower,ratio,.5_dp*dt, &
        viscosity,conduction,diffusion,barodiffusion,candidate,temps,child,child_t,transport_boundary,stage_theta,success)
      if (.not.success) return
      transport_sum=transport_sum+transport_boundary;theta=min(theta,stage_theta)
    end subroutine transport_half

    subroutine reaction_half(success)
      logical, intent(out) :: success
      integer :: c
      real(dp) :: after
      success=.false.
      do c=1,nvar
        before(c)=coarse%dx*coarse%dy*coarse%dz*sum(candidate(c,:,:,:)*coarse%volume_fraction,mask=active_coarse)+ &
          fine%dx*fine%dy*fine%dz*sum(child(c,:,:,:)*fine%volume_fraction,mask=active_fine)
      end do
      call advance_reactive_chemistry_3d(species,reactions,candidate,temps,n(1),n(2),n(3), &
        .5_dp*dt,rtol,atol,success,active_coarse,chemistry_integrator)
      if (.not.success) return
      call advance_reactive_chemistry_3d(species,reactions,child,child_t,nf(1),nf(2),nf(3), &
        .5_dp*dt,rtol,atol,success,active_fine,chemistry_integrator)
      if (.not.success) return
      do c=1,nvar
        after=coarse%dx*coarse%dy*coarse%dz*sum(candidate(c,:,:,:)*coarse%volume_fraction,mask=active_coarse)+ &
          fine%dx*fine%dy*fine%dz*sum(child(c,:,:,:)*fine%volume_fraction,mask=active_fine)
        chemical(c)=chemical(c)+after-before(c)
      end do
      call restrict_reactive_eb_state_3d(species,fine,ratio,child,child_t, &
        candidate(:,lower(1):upper(1),lower(2):upper(2),lower(3):upper(3)), &
        temps(lower(1):upper(1),lower(2):upper(2),lower(3):upper(3)),success)
    end subroutine reaction_half
  end subroutine advance_two_level_split

  ! One static fine patch, a common Euler interval and physical outflow faces.
  ! All four state arrays and both per-step inventories form one transaction.
  subroutine advance_reactive_eb_two_level_euler_3d(species,coarse,fine,lower,ratio,solver,dt, &
      coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall_change,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: coarse,fine
    integer, intent(in) :: lower(3),ratio
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:)
    real(dp), intent(inout) :: fine_state(:,:,:,:),fine_temperature(:,:,:),boundary_change(:),wall_change(:)
    logical, intent(out) :: ok
    call advance_two_level_euler(species,coarse,fine,lower,ratio,solver,dt, &
      coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall_change,ok)
  end subroutine advance_reactive_eb_two_level_euler_3d

  subroutine advance_reactive_eb_two_level_transport_3d(species,transport,coarse,fine,lower,ratio,dt, &
      viscosity,conduction,diffusion,barodiffusion,coarse_state,coarse_temperature,fine_state,fine_temperature, &
      boundary_change,minimum_theta,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(eb_geometry_3d), intent(in) :: coarse,fine
    integer, intent(in) :: lower(3),ratio
    real(dp), intent(in) :: dt
    logical, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:),fine_state(:,:,:,:),fine_temperature(:,:,:)
    real(dp), intent(inout) :: boundary_change(:),minimum_theta
    logical, intent(out) :: ok
    real(dp) :: wall(size(boundary_change))
    call advance_two_level_euler(species,coarse,fine,lower,ratio,'rusanov',dt, &
      coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall,ok, &
      transport,viscosity,conduction,diffusion,barodiffusion,minimum_theta)
  end subroutine advance_reactive_eb_two_level_transport_3d

  subroutine advance_two_level_euler(species,coarse,fine,lower,ratio,solver,dt, &
      coarse_state,coarse_temperature,fine_state,fine_temperature,boundary_change,wall_change,ok, &
      transport,viscosity,conduction,diffusion,barodiffusion,minimum_theta)
    type(nasa7_species), intent(in) :: species(:)
    type(eb_geometry_3d), intent(in) :: coarse,fine
    integer, intent(in) :: lower(3),ratio
    character(len=*), intent(in) :: solver
    real(dp), intent(in) :: dt
    real(dp), intent(inout) :: coarse_state(:,:,:,:),coarse_temperature(:,:,:)
    real(dp), intent(inout) :: fine_state(:,:,:,:),fine_temperature(:,:,:),boundary_change(:),wall_change(:)
    logical, intent(out) :: ok
    type(gas_transport_species), optional, intent(in) :: transport(:)
    logical, optional, intent(in) :: viscosity,conduction,diffusion,barodiffusion
    real(dp), optional, intent(inout) :: minimum_theta
    real(dp) :: coarse_theta,fine_theta
    logical :: faces(6)
    type(eb_geometry_3d) :: restricted,section
    type(eb_integrated_flux_3d) :: coarse_flux,fine_flux,patch_flux,effective
    real(dp), allocatable :: initial(:,:,:,:),initial_t(:,:,:),candidate(:,:,:,:),temps(:,:,:)
    real(dp), allocatable :: child(:,:,:,:),child_t(:,:,:),ghost(:,:,:,:),ghost_t(:,:,:)
    real(dp), allocatable :: coarse_wall(:,:,:,:),fine_wall(:,:,:,:),boundary(:),wall(:)
    logical, allocatable :: refined(:,:,:)
    integer :: nvar,n(3),nf(3),upper(3),c
    real(dp) :: volume
    logical :: valid
    ok=.false.
    if (.not.ieee_is_finite(dt)) return
    if (dt<=0._dp.or.ratio<1) return
    if (.not.coarse%is_valid().or..not.fine%is_valid()) return
    n=[coarse%nx,coarse%ny,coarse%nz];nf=[fine%nx,fine%ny,fine%nz]
    if (any(mod(nf,ratio)/=0).or.any(lower<1).or.any(lower>n)) return
    if (any(nf/ratio>n-lower+1)) return
    upper=lower+nf/ratio-1;nvar=reactive_nvar(size(species))
    if (any(shape(coarse_state)/=[nvar,n]).or.any(shape(coarse_temperature)/=n).or. &
      any(shape(fine_state)/=[nvar,nf]).or.any(shape(fine_temperature)/=nf)) return
    if (size(boundary_change)/=nvar.or.size(wall_change)/=nvar) return
    call coarsen_eb_geometry_3d(fine,ratio,restricted,valid)
    if (.not.valid) return
    call extract_eb_geometry_3d(coarse,lower,upper,section,valid)
    if (.not.valid) return
    if (.not.matching_metrics(section,restricted)) return
    allocate(initial,source=coarse_state);allocate(initial_t,source=coarse_temperature)
    call restrict_reactive_eb_state_3d(species,fine,ratio,fine_state,fine_temperature, &
      initial(:,lower(1):upper(1),lower(2):upper(2),lower(3):upper(3)), &
      initial_t(lower(1):upper(1),lower(2):upper(2),lower(3):upper(3)),valid)
    if (.not.valid) return
    allocate(ghost(nvar,0:nf(1)+1,0:nf(2)+1,0:nf(3)+1),ghost_t(0:nf(1)+1,0:nf(2)+1,0:nf(3)+1))
    call fill_reactive_eb_patch_halo_3d(species,coarse,fine,lower,ratio,initial,initial_t, &
      fine_state,fine_temperature,ghost,ghost_t,valid)
    if (.not.valid) return
    allocate(candidate,source=initial);allocate(temps,source=initial_t)
    allocate(child,source=fine_state);allocate(child_t,source=fine_temperature)
    if (present(transport)) then
      faces=[lower(1)>1,upper(1)<n(1),lower(2)>1,upper(2)<n(2),lower(3)>1,upper(3)<n(3)]
      call advance_reactive_eb_transport_patch_euler_3d(species,transport,coarse,dt,viscosity,conduction, &
        diffusion,barodiffusion,candidate,temps,coarse_flux,coarse_theta,valid)
      if (.not.valid) return
      call advance_reactive_eb_transport_patch_euler_3d(species,transport,fine,dt,viscosity,conduction, &
        diffusion,barodiffusion,child,child_t,fine_flux,fine_theta,valid,ghost,ghost_t,faces)
      if (.not.valid) return
    else
    call advance_reactive_eb_euler_3d(species,initial,initial_t,coarse,solver,dt,candidate,temps,valid, &
      integrated_flux=coarse_flux)
    if (.not.valid) return
    call advance_reactive_eb_euler_3d(species,fine_state,fine_temperature,fine,solver,dt,child,child_t,valid, &
      ghost,ghost_t,fine_flux)
    if (.not.valid) return
    end if
    call restrict_eb_integrated_flux_3d(fine,ratio,fine_flux,patch_flux,valid)
    if (.not.valid) return
    effective=coarse_flux
    effective%x(:,lower(1)-1:upper(1),lower(2):upper(2),lower(3):upper(3))=patch_flux%x
    effective%y(:,lower(1):upper(1),lower(2)-1:upper(2),lower(3):upper(3))=patch_flux%y
    effective%z(:,lower(1):upper(1),lower(2):upper(2),lower(3)-1:upper(3))=patch_flux%z
    allocate(refined(n(1),n(2),n(3)));refined=.false.
    refined(lower(1):upper(1),lower(2):upper(2),lower(3):upper(3))=.true.
    call reflux_reactive_eb_state_3d(species,coarse,refined,coarse_flux,effective,candidate,temps,valid)
    if (.not.valid) return
    call restrict_reactive_eb_state_3d(species,fine,ratio,child,child_t, &
      candidate(:,lower(1):upper(1),lower(2):upper(2),lower(3):upper(3)), &
      temps(lower(1):upper(1),lower(2):upper(2),lower(3):upper(3)),valid)
    if (.not.valid) return
    allocate(coarse_wall(nvar,n(1),n(2),n(3)),fine_wall(nvar,nf(1),nf(2),nf(3)),boundary(nvar),wall(nvar))
    coarse_wall=0._dp;fine_wall=0._dp
    if (.not.present(transport)) then
    call reactive_eb_slip_wall_source_3d(species,initial,initial_t,coarse,coarse_wall,valid)
    if (.not.valid) return
    call reactive_eb_slip_wall_source_3d(species,fine_state,fine_temperature,fine,fine_wall,valid)
    if (.not.valid) return
    end if
    volume=coarse%dx*coarse%dy*coarse%dz
    do c=1,nvar
      boundary(c)=coarse%dy*coarse%dz*sum( &
        coarse%x_face_fraction(0,:,:)*effective%x(c,0,:,:)-coarse%x_face_fraction(n(1),:,:)*effective%x(c,n(1),:,:)) &
        +coarse%dx*coarse%dz*sum( &
        coarse%y_face_fraction(:,0,:)*effective%y(c,:,0,:)-coarse%y_face_fraction(:,n(2),:)*effective%y(c,:,n(2),:)) &
        +coarse%dx*coarse%dy*sum( &
        coarse%z_face_fraction(:,:,0)*effective%z(c,:,:,0)-coarse%z_face_fraction(:,:,n(3))*effective%z(c,:,:,n(3)))
      wall(c)=dt*(volume*sum(coarse_wall(c,:,:,:)*coarse%volume_fraction,mask=.not.refined) &
        +fine%dx*fine%dy*fine%dz*sum(fine_wall(c,:,:,:)*fine%volume_fraction))
    end do
    if (.not.all(ieee_is_finite(boundary)).or..not.all(ieee_is_finite(wall))) return
    coarse_state=candidate;coarse_temperature=temps;fine_state=child;fine_temperature=child_t
    boundary_change=boundary;wall_change=wall
    if (present(minimum_theta)) minimum_theta=min(coarse_theta,fine_theta)
    ok=.true.
  end subroutine advance_two_level_euler

  logical function matching_metrics(a,b) result(valid)
    type(eb_geometry_3d), intent(in) :: a,b
    real(dp), parameter :: tol=2048._dp*epsilon(1._dp)
    real(dp) :: length_scale,area_scale
    valid=.false.
    if (any([a%nx,a%ny,a%nz]/=[b%nx,b%ny,b%nz])) return
    length_scale=max(1._dp,abs(a%x_lower),abs(a%x_upper),abs(a%y_lower),abs(a%y_upper), &
      abs(a%z_lower),abs(a%z_upper))
    area_scale=max(a%dx*a%dy,a%dx*a%dz,a%dy*a%dz)
    if (abs(a%x_lower-b%x_lower)>tol*length_scale.or.abs(a%x_upper-b%x_upper)>tol*length_scale) return
    if (abs(a%dx-b%dx)>tol*a%dx) return
    if (abs(a%y_lower-b%y_lower)>tol*length_scale.or.abs(a%y_upper-b%y_upper)>tol*length_scale) return
    if (abs(a%dy-b%dy)>tol*a%dy) return
    if (abs(a%z_lower-b%z_lower)>tol*length_scale.or.abs(a%z_upper-b%z_upper)>tol*length_scale) return
    if (abs(a%dz-b%dz)>tol*a%dz) return
    if (any(abs(a%volume_fraction-b%volume_fraction)>tol*1._dp)) return
    if (any(abs(a%cell_centroid_x-b%cell_centroid_x)>tol*1._dp)) return
    if (any(abs(a%cell_centroid_y-b%cell_centroid_y)>tol*1._dp)) return
    if (any(abs(a%cell_centroid_z-b%cell_centroid_z)>tol*1._dp)) return
    if (any(abs(a%x_face_fraction-b%x_face_fraction)>tol*1._dp)) return
    if (any(abs(a%y_face_fraction-b%y_face_fraction)>tol*1._dp)) return
    if (any(abs(a%z_face_fraction-b%z_face_fraction)>tol*1._dp)) return
    if (any(abs(a%x_face_centroid_y-b%x_face_centroid_y)>tol*1._dp)) return
    if (any(abs(a%x_face_centroid_z-b%x_face_centroid_z)>tol*1._dp)) return
    if (any(abs(a%y_face_centroid_x-b%y_face_centroid_x)>tol*1._dp)) return
    if (any(abs(a%y_face_centroid_z-b%y_face_centroid_z)>tol*1._dp)) return
    if (any(abs(a%z_face_centroid_x-b%z_face_centroid_x)>tol*1._dp)) return
    if (any(abs(a%z_face_centroid_y-b%z_face_centroid_y)>tol*1._dp)) return
    if (any(abs(a%boundary_area-b%boundary_area)>tol*area_scale)) return
    if (any(abs(a%boundary_centroid_x-b%boundary_centroid_x)>tol*length_scale)) return
    if (any(abs(a%boundary_centroid_y-b%boundary_centroid_y)>tol*length_scale)) return
    if (any(abs(a%boundary_centroid_z-b%boundary_centroid_z)>tol*length_scale)) return
    if (any(abs(a%boundary_normal_x-b%boundary_normal_x)>tol*1._dp)) return
    if (any(abs(a%boundary_normal_y-b%boundary_normal_y)>tol*1._dp)) return
    if (any(abs(a%boundary_normal_z-b%boundary_normal_z)>tol*1._dp)) return
    if (any(abs(a%boundary_normal_integral_x-b%boundary_normal_integral_x)>tol*area_scale)) return
    if (any(abs(a%boundary_normal_integral_y-b%boundary_normal_integral_y)>tol*area_scale)) return
    if (any(abs(a%boundary_normal_integral_z-b%boundary_normal_integral_z)>tol*area_scale)) return
    if (any(a%cell_type/=b%cell_type)) return
    valid=.true.
  end function matching_metrics
end module eb_reactive_two_level_3d_mod
