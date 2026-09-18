module les_smagorinsky_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use simulation_config_reactive_3d_mod, only: fill_reactive_boundary_ghosts_3d, &
    valid_reactive_boundaries_3d,valid_wall_velocity_3d,fill_pressure_outlet_ghosts_3d
  use nasa7_thermo_mod, only: nasa7_species, nasa7_mass_properties
  use mixture_thermo_mod, only: mixture_mass_properties
  use reactive_1d_mod, only: reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  implicit none
  private
  public :: wale_stress,vreman_stress
  public :: compute_les_ghosted_rhs_3d
  public :: smagorinsky_stress, compute_les_rhs_3d, advance_les_ssprk2_3d
contains
  subroutine smagorinsky_stress(gradient, density, spacing, cs, viscosity, stress, ok,ci)
    real(dp), intent(in) :: gradient(3,3), density, spacing(3), cs
    real(dp), intent(out) :: viscosity, stress(3,3)
    real(dp), optional, intent(in) :: ci
    logical, intent(out) :: ok
    real(dp) :: strain(3,3), delta, trace
    integer :: d
    viscosity=0.0_dp; stress=0.0_dp; ok=.false.
    if (.not. all(ieee_is_finite(gradient))) return
    if (.not. all(ieee_is_finite([density,spacing,cs]))) return
    if (maxval(abs(gradient)) > 1.0e30_dp) return
    if (density <= 0.0_dp .or. density > 1.0e8_dp) return
    if (any(spacing < 1.0e-12_dp) .or. any(spacing > 1.0e9_dp)) return
    if (cs < 0.0_dp .or. cs > 2.0_dp) return
    delta=product(spacing)**(1.0_dp/3.0_dp)
    strain=0.5_dp*(gradient+transpose(gradient))
    trace=(strain(1,1)+strain(2,2)+strain(3,3))/3.0_dp
    viscosity=density*(cs*delta)**2*sqrt(2.0_dp*sum(strain**2))
    do d=1,3
      strain(d,d)=strain(d,d)-trace
    end do
    ! Isotropic stress is added separately, without changing eddy viscosity.
    stress=2.0_dp*viscosity*strain
    ok=.true.
    if (ok) call add_isotropic_stress(gradient,density,spacing,ci,viscosity,stress,ok)
  end subroutine

  ! Pinned PeleC LES.H WALE/Vreman closures; CI defaults to zero.
  ! Gradient rows are velocity components, columns are spatial derivatives.
  subroutine wale_stress(gradient,density,spacing,cw,viscosity,stress,ok,ci)
    real(dp), intent(in) :: gradient(3,3),density,spacing(3),cw
    real(dp), intent(out) :: viscosity,stress(3,3)
    real(dp), optional, intent(in) :: ci
    logical, intent(out) :: ok
    call algebraic_stress(gradient,density,spacing,cw,.true.,viscosity,stress,ok)
    if (ok) call add_isotropic_stress(gradient,density,spacing,ci,viscosity,stress,ok)
  end subroutine

  subroutine vreman_stress(gradient,density,spacing,cs,viscosity,stress,ok,ci)
    real(dp), intent(in) :: gradient(3,3),density,spacing(3),cs
    real(dp), intent(out) :: viscosity,stress(3,3)
    real(dp), optional, intent(in) :: ci
    logical, intent(out) :: ok
    call algebraic_stress(gradient,density,spacing,cs,.false.,viscosity,stress,ok)
    if (ok) call add_isotropic_stress(gradient,density,spacing,ci,viscosity,stress,ok)
  end subroutine

  ! All three pinned static closures use alpha=2*rho*delta**2*|S|**2.
  subroutine add_isotropic_stress(gradient,density,spacing,ci,viscosity,stress,ok)
    real(dp), intent(in) :: gradient(3,3),density,spacing(3)
    real(dp), optional, intent(in) :: ci
    real(dp), intent(inout) :: viscosity,stress(3,3)
    logical, intent(out) :: ok
    real(dp) :: pressure
    integer :: d
    ok=.true.
    if (.not.present(ci)) return
    ok=.false.
    if (.not.ieee_is_finite(ci)) then
      viscosity=0;stress=0;return
    end if
    if (ci<0.or.ci>1) then
      viscosity=0;stress=0;return
    end if
    ok=.true.
    if (ci==0) return
    pressure=(4._dp/3._dp)*ci*density*product(spacing)**(2._dp/3._dp)* &
      sum((.5_dp*(gradient+transpose(gradient)))**2)
    do d=1,3
      stress(d,d)=stress(d,d)-pressure
    end do
  end subroutine

  subroutine algebraic_stress(gradient,density,spacing,coefficient,wale,viscosity,stress,ok)
    real(dp), intent(in) :: gradient(3,3),density,spacing(3),coefficient
    logical, intent(in) :: wale
    real(dp), intent(out) :: viscosity,stress(3,3)
    logical, intent(out) :: ok
    real(dp) :: g(3,3),strain(3,3),square(3,3),deviator(3,3),scale,delta,trace,s2,d2,rate,beta
    integer :: i,j,k,l
    viscosity=0;stress=0;ok=.false.
    if (.not.all(ieee_is_finite(gradient))) return
    if (.not.all(ieee_is_finite([density,spacing,coefficient]))) return
    scale=maxval(abs(gradient))
    if (scale>1e30_dp.or.density<=0.or.density>1e8_dp) return
    if (any(spacing<1e-12_dp).or.any(spacing>1e9_dp)) return
    if (coefficient<0.or.coefficient>2) return
    ok=.true.
    if (scale==0.or.coefficient==0) return
    ! Homogeneity permits normalized invariants without fifth-power overflow.
    g=gradient/scale;strain=.5_dp*(g+transpose(g))
    delta=product(spacing)**(1._dp/3._dp)
    if (wale) then
      square=matmul(g,g);deviator=.5_dp*(square+transpose(square))
      trace=(square(1,1)+square(2,2)+square(3,3))/3
      do i=1,3
        deviator(i,i)=deviator(i,i)-trace
      end do
      s2=sum(strain**2);d2=sum(deviator**2)
      rate=0
      if (d2>0) rate=d2**1.5_dp/(s2**2.5_dp+d2**1.25_dp)
    else
      ! Sum squared 2x2 minors equals the Gram-matrix invariant B_beta;
      ! this form avoids cancellation to a negative value for rank-one shear.
      beta=0
      do i=1,2
        do j=i+1,3
          do k=1,2
            do l=k+1,3
              beta=beta+(g(k,i)*g(l,j)-g(l,i)*g(k,j))**2
            end do
          end do
        end do
      end do
      rate=2.5_dp*sqrt(beta/sum(g**2))
    end if
    viscosity=density*(coefficient*delta)**2*scale*rate
    trace=(strain(1,1)+strain(2,2)+strain(3,3))/3
    do i=1,3
      strain(i,i)=strain(i,i)-trace
    end do
    stress=2*viscosity*scale*strain
  end subroutine

  subroutine compute_les_rhs_3d(species,state,temperature,spacing,cs,prt,sct,rhs,maximum_diffusivity,ok, &
      boundary_faces,boundary_rate,inflow_state,inflow_temperature,wall_temperature,wall_velocity,model,ci,outlet_pressure)
    real(dp), optional, intent(in) :: outlet_pressure(6)
    real(dp), optional, intent(in) :: ci
    character(len=*), optional, intent(in) :: model
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),spacing(3),cs,prt,sct
    real(dp), intent(out) :: rhs(:,:,:,:),maximum_diffusivity
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    character(len=*), intent(in), optional :: boundary_faces(6)
    real(dp), intent(out), optional :: boundary_rate(:)
    character(len=32) :: faces(6)
    real(dp), allocatable :: ghost(:,:,:,:),gt(:,:,:),grhs(:,:,:,:)
    real(dp) :: rate(size(state,1))
    integer :: n(3),dims(3),f
    ok=.false.; rhs=0; maximum_diffusivity=0; faces='periodic'
    if (present(boundary_rate)) then
      boundary_rate=0
      if (size(boundary_rate)/=size(state,1)) return
    end if
    if (present(boundary_faces)) faces=boundary_faces
    if (.not.valid_reactive_boundaries_3d(faces)) return
    if (present(outlet_pressure)) then
      if (.not.all(ieee_is_finite(outlet_pressure))) return
      if (any(outlet_pressure<0).or.any(outlet_pressure>1e15_dp)) return
      if (any(outlet_pressure>0.and.faces/='outflow')) return
    end if
    do f=1,6
      select case(trim(faces(f)))
      case('periodic','outflow','inflow','slip_wall','no_slip_wall','isothermal_wall')
      case default
        return
      end select
    end do
    if (present(wall_velocity)) then
      if (.not.valid_wall_velocity_3d(faces,wall_velocity)) return
    end if
    if (any(faces=='isothermal_wall')) then
      if (.not.present(wall_temperature)) return
      if (size(wall_temperature)/=6) return
      do f=1,6
        if (faces(f)/='isothermal_wall') cycle
        if (.not.ieee_is_finite(wall_temperature(f))) return
        if (wall_temperature(f)<=0) return
      end do
    end if
    if (all(faces=='periodic')) then
      call compute_les_face_rhs_3d(species,state,temperature,spacing,cs,prt,sct,rhs,maximum_diffusivity,ok,model=model,ci=ci)
      return
    end if
    n=shape(temperature); dims=n+4
    if (any(shape(rhs)/=shape(state)).or.any(n<2)) return
    allocate(ghost(size(state,1),dims(1),dims(2),dims(3)),gt(dims(1),dims(2),dims(3)))
    allocate(grhs,mold=ghost)
    call fill_reactive_boundary_ghosts_3d(state,temperature,faces,ghost,gt,ok,viscous_wall=.true., &
      inflow_state=inflow_state,inflow_temperature=inflow_temperature,wall_velocity=wall_velocity)
    if (.not.ok) return
    call fill_pressure_outlet_ghosts_3d(species,n,faces,outlet_pressure,ghost,gt,ok)
    if (.not.ok) return
    call compute_les_face_rhs_3d(species,ghost,gt,spacing,cs,prt,sct,grhs,maximum_diffusivity,ok,n,rate, &
      faces,wall_temperature,model=model,ci=ci)
    if (.not.ok) return
    rhs=grhs(:,3:n(1)+2,3:n(2)+2,3:n(3)+2)
    if (present(boundary_rate)) boundary_rate=rate
  end subroutine

  ! Two supplied halo layers on every face; return only the active-cell RHS.
  subroutine compute_les_ghosted_rhs_3d(species,state,temperature,cells,spacing,cs,prt,sct,rhs,diffusivity,ok, &
      boundary_faces,boundary_rate,wall_temperature,boundary_face_mask,model,ci,fx,fy,fz)
    ! Flux density in the positive coordinate direction; all three or none.
    ! Face arrays publish together only after a successful RHS evaluation.
    real(dp), optional, intent(inout) :: fx(:,0:,:,:),fy(:,:,0:,:),fz(:,:,:,0:)
    real(dp), allocatable :: xflux(:,:,:,:),yflux(:,:,:,:),zflux(:,:,:,:)
    real(dp), optional, intent(in) :: ci
    character(len=*), optional, intent(in) :: model
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),spacing(3),cs,prt,sct
    integer, intent(in) :: cells(3)
    real(dp), intent(out) :: rhs(:,:,:,:),diffusivity
    character(len=*), optional, intent(in) :: boundary_faces(6)
    real(dp), optional, intent(in) :: wall_temperature(:)
    real(dp), optional, intent(out) :: boundary_rate(:)
    logical, optional, intent(in) :: boundary_face_mask(6)
    integer :: f
    logical, intent(out) :: ok
    real(dp), allocatable :: full_rhs(:,:,:,:)
    ok=.false.;rhs=0;diffusivity=0
    if (present(boundary_rate)) then
      boundary_rate=0
      if (size(boundary_rate)/=size(state,1)) return
    end if
    if (present(boundary_faces)) then
      if (any(boundary_faces=='isothermal_wall')) then
        if (.not.present(wall_temperature)) return
        if (size(wall_temperature)/=6) return
        do f=1,6
          if (boundary_faces(f)/='isothermal_wall') cycle
          if (.not.ieee_is_finite(wall_temperature(f))) return
          if (wall_temperature(f)<=0) return
        end do
      end if
    end if
    if (any(cells<2).or.any(cells>huge(1)-4)) return
    if (any(shape(temperature)/=cells+4).or.any(shape(rhs)/=[size(state,1),cells])) return
    if ((present(fx).neqv.present(fy)).or.(present(fx).neqv.present(fz))) return
    if (present(fx)) then
      if (any(shape(fx)/=[size(state,1),cells(1)+1,cells(2),cells(3)])) return
      if (any(shape(fy)/=[size(state,1),cells(1),cells(2)+1,cells(3)])) return
      if (any(shape(fz)/=[size(state,1),cells(1),cells(2),cells(3)+1])) return
      allocate(xflux(size(state,1),0:cells(1),cells(2),cells(3)))
      allocate(yflux(size(state,1),cells(1),0:cells(2),cells(3)))
      allocate(zflux(size(state,1),cells(1),cells(2),0:cells(3)))
    end if
    allocate(full_rhs,mold=state)
    call compute_les_face_rhs_3d(species,state,temperature,spacing,cs,prt,sct,full_rhs,diffusivity,ok,cells, &
      boundary_rate,boundary_faces,wall_temperature,boundary_face_mask,model,ci,xflux,yflux,zflux)
    if (ok) then
      rhs=full_rhs(:,3:cells(1)+2,3:cells(2)+2,3:cells(3)+2)
      if (present(fx)) then
        fx=xflux;fy=yflux;fz=zflux
      end if
    end if
  end subroutine

  subroutine compute_les_face_rhs_3d(species,state,temperature,spacing,cs,prt,sct,rhs,maximum_diffusivity,ok, &
      active_cells,boundary_rate,boundary_faces,wall_temperature,boundary_face_mask,model,ci,fx,fy,fz)
    real(dp), optional, intent(out) :: fx(:,0:,:,:),fy(:,:,0:,:),fz(:,:,:,0:)
    real(dp), optional, intent(in) :: ci
    character(len=*), optional, intent(in) :: model
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),spacing(3),cs,prt,sct
    real(dp), intent(out) :: rhs(:,:,:,:), maximum_diffusivity
    logical, intent(out) :: ok
    integer, intent(in), optional :: active_cells(3)
    character(len=*), intent(in), optional :: boundary_faces(6)
    real(dp), intent(in), optional :: wall_temperature(:)
    logical, optional, intent(in) :: boundary_face_mask(6)
    logical :: include_face(6)
    real(dp), intent(out), optional :: boundary_rate(:)
    real(dp) :: rate(size(state,1))
    character(len=16) :: selected_model
    logical :: active_face
    real(dp), allocatable :: primitive(:,:,:,:), temps(:,:,:), cp(:,:,:), cv(:,:,:), enthalpies(:,:,:,:)
    real(dp) :: q(reactive_nprim(size(species))), y(size(species)), grad(3,3), tau(3,3), flux(size(state,1))
    real(dp) :: jflux(size(species)),hav(size(species)),r,mw,gamma,h,e,s,cpg,cvg,sound,t,mu,rho,heat_flux
    real(dp) :: face_y(size(species)),uface(3),cpface,cvface,diffusivity,cellvolume,isotropic,iso_diff
    integer :: dims(3),l(3),ridx(3),lo(3),hi(3),lo2(3),hi2(3),i,j,k,n,a,b,wall_face
    logical :: local_ok
    rhs=0.0_dp; maximum_diffusivity=0.0_dp; ok=.false.; rate=0
    isotropic=0
    if (present(ci)) isotropic=ci
    if (.not.ieee_is_finite(isotropic)) return
    if (isotropic<0.or.isotropic>1) return
    if (present(fx)) then
      fx=0;fy=0;fz=0
    end if
    selected_model='smagorinsky'
    if (present(model)) then
      select case(model)
      case('smagorinsky','wale','vreman')
        selected_model=model
      case default
        return
      end select
    end if
    include_face=.true.
    if (present(boundary_face_mask)) include_face=boundary_face_mask
    if (present(boundary_rate)) boundary_rate=0
    dims=shape(temperature)
    if (any(dims < 2)) return
    if (size(state,1) /= reactive_nvar(size(species)) .or. size(species) < 2) return
    if (any(shape(state(1,:,:,:)) /= dims) .or. any(shape(rhs) /= shape(state))) return
    if (.not. all(ieee_is_finite([spacing,cs,prt,sct]))) return
    if (any(spacing < 1.0e-12_dp) .or. any(spacing > 1.0e9_dp)) return
    if (cs < 0.0_dp .or. cs > 2.0_dp .or. min(prt,sct) < 1.0e-6_dp .or. max(prt,sct) > 1.0e6_dp) return
    allocate(primitive(size(q),dims(1),dims(2),dims(3)),temps(dims(1),dims(2),dims(3)))
    allocate(cp(dims(1),dims(2),dims(3)),cv(dims(1),dims(2),dims(3)))
    allocate(enthalpies(size(species),dims(1),dims(2),dims(3)))
    do k=1,dims(3)
      do j=1,dims(2)
        do i=1,dims(1)
          call reactive_conserved_to_primitive(species,state(:,i,j,k),temperature(i,j,k),q,t,sound,local_ok)
          if (.not. local_ok) return
          if (maxval(abs(q(2:4))) > 1.0e7_dp .or. q(1) > 1.0e8_dp) return
          primitive(:,i,j,k)=q; temps(i,j,k)=t; y=q(6:)
          call mixture_mass_properties(species,y,t,mw,r,cpg,cvg,gamma,h,e,s,local_ok)
          if (.not. local_ok) return
          cp(i,j,k)=cpg; cv(i,j,k)=cvg
          do n=1,size(species)
            call nasa7_mass_properties(species(n),t,cpg,cvg,h,e,s,local_ok)
            if (.not. local_ok) return
            enthalpies(n,i,j,k)=h
          end do
        end do
      end do
    end do
    if (cs == 0.0_dp.and.isotropic==0) then
      ok=.true.; return
    end if
    cellvolume=product(spacing)
    do a=1,3
      do k=1,dims(3)
        do j=1,dims(2)
          do i=1,dims(1)
            l=[i,j,k]
            if (present(active_cells)) then
              active_face=l(a)>=2.and.l(a)<=active_cells(a)+2
              do b=1,3
                if (b==a) cycle
                active_face=active_face.and.l(b)>=3.and.l(b)<=active_cells(b)+2
              end do
              if (.not.active_face) cycle
            end if
            ridx=l; ridx(a)=modulo(l(a),dims(a))+1
            grad(:,a)=(primitive(2:4,ridx(1),ridx(2),ridx(3))-primitive(2:4,i,j,k))/spacing(a)
            do b=1,3
              if (b == a) cycle
              lo=l; hi=l; lo2=ridx; hi2=ridx
              lo(b)=modulo(l(b)-2,dims(b))+1; hi(b)=modulo(l(b),dims(b))+1
              lo2(b)=modulo(ridx(b)-2,dims(b))+1; hi2(b)=modulo(ridx(b),dims(b))+1
              grad(:,b)=(primitive(2:4,hi(1),hi(2),hi(3))-primitive(2:4,lo(1),lo(2),lo(3)) + &
                primitive(2:4,hi2(1),hi2(2),hi2(3))-primitive(2:4,lo2(1),lo2(2),lo2(3)))/(4.0_dp*spacing(b))
            end do
            rho=0.5_dp*(primitive(1,i,j,k)+primitive(1,ridx(1),ridx(2),ridx(3)))
            select case(selected_model)
            case('smagorinsky')
              call smagorinsky_stress(grad,rho,spacing,cs,mu,tau,local_ok,ci)
            case('wale')
              call wale_stress(grad,rho,spacing,cs,mu,tau,local_ok,ci)
            case('vreman')
              call vreman_stress(grad,rho,spacing,cs,mu,tau,local_ok,ci)
            end select
            if (.not. local_ok) then
              rhs=0.0_dp; maximum_diffusivity=0.0_dp; return
            end if
            cpface=0.5_dp*(cp(i,j,k)+cp(ridx(1),ridx(2),ridx(3)))
            cvface=0.5_dp*(cv(i,j,k)+cv(ridx(1),ridx(2),ridx(3)))
            uface=0.5_dp*(primitive(2:4,i,j,k)+primitive(2:4,ridx(1),ridx(2),ridx(3)))
            face_y=0.5_dp*(primitive(6:,i,j,k)+primitive(6:,ridx(1),ridx(2),ridx(3)))
            hav=0.5_dp*(enthalpies(:,i,j,k)+enthalpies(:,ridx(1),ridx(2),ridx(3)))
            jflux=-mu/sct*(primitive(6:,ridx(1),ridx(2),ridx(3))-primitive(6:,i,j,k))/spacing(a)
            jflux=jflux-face_y*sum(jflux)
            heat_flux=-mu*cpface/prt*(temps(ridx(1),ridx(2),ridx(3))-temps(i,j,k))/spacing(a)
            wall_face=0
            if (present(active_cells).and.present(boundary_faces).and.present(wall_temperature)) then
              if (l(a)==2) wall_face=2*a-1
              if (l(a)==active_cells(a)+2) wall_face=2*a
              if (wall_face>0) then
                if (boundary_faces(wall_face)=='isothermal_wall') then
                  ! Wall-to-cell distance is half a cell; ghost energy stays admissible.
                  if (modulo(wall_face,2)==1) then
                    heat_flux=-2*mu*cpface/prt*(temps(ridx(1),ridx(2),ridx(3))-wall_temperature(wall_face))/spacing(a)
                  else
                    heat_flux=-2*mu*cpface/prt*(wall_temperature(wall_face)-temps(i,j,k))/spacing(a)
                  end if
                else
                  wall_face=0
                end if
              end if
            end if
            flux=0.0_dp; flux(2:4)=-tau(:,a)
            flux(5)=-dot_product(uface,tau(:,a))+heat_flux+dot_product(hav,jflux)
            flux(6:)=jflux
            if (present(fx)) then
              select case(a)
              case(1)
                fx(:,i-2,j-2,k-2)=flux
              case(2)
                fy(:,i-2,j-2,k-2)=flux
              case(3)
                fz(:,i-2,j-2,k-2)=flux
              end select
            end if
            if (present(active_cells)) then
              if (l(a)==2.and.include_face(2*a-1)) rate=rate+flux*cellvolume/spacing(a)
              if (l(a)==active_cells(a)+2.and.include_face(2*a)) rate=rate-flux*cellvolume/spacing(a)
            end if
            ! One face flux, equal and opposite for both cells: exact telescoping.
            rhs(:,i,j,k)=rhs(:,i,j,k)-flux/spacing(a)
            rhs(:,ridx(1),ridx(2),ridx(3))=rhs(:,ridx(1),ridx(2),ridx(3))+flux/spacing(a)
            diffusivity=mu*max(4.0_dp/3.0_dp,cpface/(cvface*prt),1.0_dp/sct) / &
              min(primitive(1,i,j,k),primitive(1,ridx(1),ridx(2),ridx(3)))
            if (isotropic>0) then
              ! Additional frozen-gradient estimate; not a nonlinear stability proof.
              iso_diff=4*isotropic*rho*product(spacing)**(2._dp/3._dp)* &
                sqrt(2*sum((.5_dp*(grad+transpose(grad)))**2))/ &
                min(primitive(1,i,j,k),primitive(1,ridx(1),ridx(2),ridx(3)))
              diffusivity=diffusivity+iso_diff
            end if
            if (wall_face>0) diffusivity=2*diffusivity
            maximum_diffusivity=max(maximum_diffusivity,diffusivity)
          end do
        end do
      end do
    end do
    if (present(boundary_rate)) boundary_rate=rate
    ok=all(ieee_is_finite(rhs)) .and. ieee_is_finite(maximum_diffusivity) .and. cellvolume > 0.0_dp
    if (.not. ok) then
      rhs=0.0_dp; maximum_diffusivity=0.0_dp
    end if
  end subroutine

  subroutine advance_les_ssprk2_3d(species,state,temperature,spacing,cs,prt,sct,dt,ok, &
      boundary_faces,boundary_change,inflow_state,inflow_temperature,wall_temperature,wall_velocity,model,ci,outlet_pressure)
    real(dp), optional, intent(in) :: outlet_pressure(6)
    real(dp), optional, intent(in) :: ci
    character(len=*), optional, intent(in) :: model
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp), intent(in) :: spacing(3),cs,prt,sct,dt
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: inflow_state(:,:),inflow_temperature(:),wall_temperature(:),wall_velocity(:,:)
    character(len=*), intent(in), optional :: boundary_faces(6)
    real(dp), intent(out), optional :: boundary_change(:)
    real(dp) :: rate1(size(state,1)),rate2(size(state,1))
    real(dp), allocatable :: rhs(:,:,:,:),stage(:,:,:,:),candidate(:,:,:,:),t1(:,:,:),t2(:,:,:)
    real(dp) :: diffusivity,isotropic
    integer :: dims(3)
    logical :: local_ok
    isotropic=0
    if (present(ci)) isotropic=ci
    ok=.false.; dims=shape(temperature)
    if (present(boundary_change)) then
      boundary_change=0
      if (size(boundary_change)/=size(state,1)) return
    end if
    if (.not. ieee_is_finite(dt)) return
    if (dt < 0.0_dp .or. dt > 1.0e6_dp) return
    allocate(rhs,mold=state); allocate(stage,mold=state); allocate(candidate,mold=state)
    allocate(t1,mold=temperature); allocate(t2,mold=temperature)
    call compute_les_rhs_3d(species,state,temperature,spacing,cs,prt,sct,rhs,diffusivity,local_ok,boundary_faces,rate1, &
      inflow_state,inflow_temperature,wall_temperature,wall_velocity,model,ci,outlet_pressure)
    if (.not. local_ok) return
    if (dt*diffusivity*sum(1.0_dp/spacing**2) > 0.2_dp) return
    if ((cs==0.and.isotropic==0).or.dt==0) then
      ok=.true.; return
    end if
    stage=state+dt*rhs
    call recover_reactive_temperatures_3d(species,stage,temperature,dims(1),dims(2),dims(3),t1,local_ok)
    if (.not. local_ok) return
    call compute_les_rhs_3d(species,stage,t1,spacing,cs,prt,sct,rhs,diffusivity,local_ok,boundary_faces,rate2, &
      inflow_state,inflow_temperature,wall_temperature,wall_velocity,model,ci,outlet_pressure)
    if (.not. local_ok) return
    if (dt*diffusivity*sum(1.0_dp/spacing**2) > 0.2_dp) return
    candidate=0.5_dp*(state+stage+dt*rhs)
    call recover_reactive_temperatures_3d(species,candidate,t1,dims(1),dims(2),dims(3),t2,local_ok)
    if (.not. local_ok) return
    if (present(boundary_change)) boundary_change=0.5_dp*dt*(rate1+rate2)
    state=candidate; temperature=t2; ok=.true.
  end subroutine
end module les_smagorinsky_3d_mod
