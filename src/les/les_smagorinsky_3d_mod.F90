module les_smagorinsky_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species, nasa7_mass_properties
  use mixture_thermo_mod, only: mixture_mass_properties
  use reactive_1d_mod, only: reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  implicit none
  private
  public :: smagorinsky_stress, compute_les_rhs_3d, advance_les_ssprk2_3d
contains
  subroutine smagorinsky_stress(gradient, density, spacing, cs, viscosity, stress, ok)
    real(dp), intent(in) :: gradient(3,3), density, spacing(3), cs
    real(dp), intent(out) :: viscosity, stress(3,3)
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
    ! Deviatoric compressible Smagorinsky only: no modeled isotropic SGS energy.
    stress=2.0_dp*viscosity*strain
    ok=.true.
  end subroutine

  subroutine compute_les_rhs_3d(species, state, temperature, spacing, cs, prt, sct, rhs, maximum_diffusivity, ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),spacing(3),cs,prt,sct
    real(dp), intent(out) :: rhs(:,:,:,:), maximum_diffusivity
    logical, intent(out) :: ok
    real(dp), allocatable :: primitive(:,:,:,:), temps(:,:,:), cp(:,:,:), cv(:,:,:), enthalpies(:,:,:,:)
    real(dp) :: q(reactive_nprim(size(species))), y(size(species)), grad(3,3), tau(3,3), flux(size(state,1))
    real(dp) :: jflux(size(species)),hav(size(species)),r,mw,gamma,h,e,s,cpg,cvg,sound,t,mu,rho,heat_flux
    real(dp) :: face_y(size(species)),uface(3),cpface,cvface,diffusivity,cellvolume
    integer :: dims(3),l(3),ridx(3),lo(3),hi(3),lo2(3),hi2(3),i,j,k,n,a,b
    logical :: local_ok
    rhs=0.0_dp; maximum_diffusivity=0.0_dp; ok=.false.
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
    if (cs == 0.0_dp) then
      ok=.true.; return
    end if
    cellvolume=product(spacing)
    do a=1,3
      do k=1,dims(3)
        do j=1,dims(2)
          do i=1,dims(1)
            l=[i,j,k]; ridx=l; ridx(a)=modulo(l(a),dims(a))+1
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
            call smagorinsky_stress(grad,rho,spacing,cs,mu,tau,local_ok)
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
            flux=0.0_dp; flux(2:4)=-tau(:,a)
            flux(5)=-dot_product(uface,tau(:,a))+heat_flux+dot_product(hav,jflux)
            flux(6:)=jflux
            ! One face flux, equal and opposite for both cells: exact telescoping.
            rhs(:,i,j,k)=rhs(:,i,j,k)-flux/spacing(a)
            rhs(:,ridx(1),ridx(2),ridx(3))=rhs(:,ridx(1),ridx(2),ridx(3))+flux/spacing(a)
            diffusivity=mu*max(4.0_dp/3.0_dp,cpface/(cvface*prt),1.0_dp/sct) / &
              min(primitive(1,i,j,k),primitive(1,ridx(1),ridx(2),ridx(3)))
            maximum_diffusivity=max(maximum_diffusivity,diffusivity)
          end do
        end do
      end do
    end do
    ok=all(ieee_is_finite(rhs)) .and. ieee_is_finite(maximum_diffusivity) .and. cellvolume > 0.0_dp
    if (.not. ok) then
      rhs=0.0_dp; maximum_diffusivity=0.0_dp
    end if
  end subroutine

  subroutine advance_les_ssprk2_3d(species,state,temperature,spacing,cs,prt,sct,dt,ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp), intent(in) :: spacing(3),cs,prt,sct,dt
    logical, intent(out) :: ok
    real(dp), allocatable :: rhs(:,:,:,:),stage(:,:,:,:),candidate(:,:,:,:),t1(:,:,:),t2(:,:,:)
    real(dp) :: diffusivity
    integer :: dims(3)
    logical :: local_ok
    ok=.false.; dims=shape(temperature)
    if (.not. ieee_is_finite(dt)) return
    if (dt < 0.0_dp .or. dt > 1.0e6_dp) return
    allocate(rhs,mold=state); allocate(stage,mold=state); allocate(candidate,mold=state)
    allocate(t1,mold=temperature); allocate(t2,mold=temperature)
    call compute_les_rhs_3d(species,state,temperature,spacing,cs,prt,sct,rhs,diffusivity,local_ok)
    if (.not. local_ok) return
    if (dt*diffusivity*sum(1.0_dp/spacing**2) > 0.2_dp) return
    if (cs == 0.0_dp .or. dt == 0.0_dp) then
      ok=.true.; return
    end if
    stage=state+dt*rhs
    call recover_reactive_temperatures_3d(species,stage,temperature,dims(1),dims(2),dims(3),t1,local_ok)
    if (.not. local_ok) return
    call compute_les_rhs_3d(species,stage,t1,spacing,cs,prt,sct,rhs,diffusivity,local_ok)
    if (.not. local_ok) return
    if (dt*diffusivity*sum(1.0_dp/spacing**2) > 0.2_dp) return
    candidate=0.5_dp*(state+stage+dt*rhs)
    call recover_reactive_temperatures_3d(species,candidate,t1,dims(1),dims(2),dims(3),t2,local_ok)
    if (.not. local_ok) return
    state=candidate; temperature=t2; ok=.true.
  end subroutine
end module les_smagorinsky_3d_mod
