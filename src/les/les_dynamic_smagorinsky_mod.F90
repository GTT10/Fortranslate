module les_dynamic_smagorinsky_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use les_smagorinsky_3d_mod, only: smagorinsky_stress
  use les_filter_3d_mod, only: les_filter_weights,filter_les_field_3d
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: mixture_mass_properties
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_conserved_to_primitive
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  implicit none
  private
  public :: les_face_mean
  public :: advance_dynamic_les_ssprk2_3d
  public :: compute_dynamic_les_rhs_3d
  public :: compute_dynamic_les_rhs_ghosted_3d
  public :: compute_dynamic_les_coefficients_3d
  public :: dynamic_smagorinsky_quantities,dynamic_smagorinsky_coefficients
  ! Pinned PeleC Constants.H small_num(), including its coefficient floor.
  real(dp), parameter :: regularizer=1e-8_dp
contains
  ! Caller supplies ng valid donor cells on every side, all at the patch spacing.
  ! The periodic workspace seam is outside the complete filter/gradient stencil.
  subroutine compute_dynamic_les_rhs_ghosted_3d(species,state,temperature,spacing,test_model,test_ratio,ng, &
      rhs,fx,fy,fz,ok,grid_ratio,harmonic_mean,maximum_diffusivity)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),spacing(3)
    character(len=*), intent(in) :: test_model
    integer, intent(in) :: test_ratio,ng
    real(dp), intent(inout) :: rhs(:,:,:,:),fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: grid_ratio
    logical, optional, intent(in) :: harmonic_mean
    real(dp), optional, intent(inout) :: maximum_diffusivity
    real(dp), allocatable :: weights(:),full_rhs(:,:,:,:),x(:,:,:,:),y(:,:,:,:),z(:,:,:,:)
    integer :: n(3),padded(3),nv,lo(3),hi(3)
    real(dp) :: diffusion
    logical :: valid
    ok=.false.;nv=reactive_nvar(size(species))
    n=[size(rhs,2),size(rhs,3),size(rhs,4)]
    if(size(species)<1.or.size(rhs,1)/=nv.or.any(n<1)) return
    call les_filter_weights(test_model,test_ratio,weights,valid)
    if(.not.valid) return
    ! Test filtering radius + one derivative + three smoothing cells + one face neighbor.
    if(ng<(size(weights)-1)/2+5.or.ng>128) return
    padded=shape(temperature)
    if(any(padded/=n+2*ng).or.any(shape(state)/=[nv,padded])) return
    if(any(shape(fx)/=[nv,n(1)+1,n(2),n(3)])) return
    if(any(shape(fy)/=[nv,n(1),n(2)+1,n(3)])) return
    if(any(shape(fz)/=[nv,n(1),n(2),n(3)+1])) return
    allocate(full_rhs(nv,padded(1),padded(2),padded(3)))
    allocate(x(nv,padded(1)+1,padded(2),padded(3)),y(nv,padded(1),padded(2)+1,padded(3)))
    allocate(z(nv,padded(1),padded(2),padded(3)+1))
    call compute_dynamic_les_rhs_3d(species,state,temperature,spacing,test_model,test_ratio,full_rhs,valid, &
      grid_ratio=grid_ratio,harmonic_mean=harmonic_mean,fx=x,fy=y,fz=z,maximum_diffusivity=diffusion)
    if(.not.valid) return
    lo=ng+1;hi=ng+n
    rhs=full_rhs(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3))
    fx=x(:,lo(1):hi(1)+1,lo(2):hi(2),lo(3):hi(3))
    fy=y(:,lo(1):hi(1),lo(2):hi(2)+1,lo(3):hi(3))
    fz=z(:,lo(1):hi(1),lo(2):hi(2),lo(3):hi(3)+1)
    ! Conservative upper bound includes the padded workspace as well as the patch.
    if(present(maximum_diffusivity)) maximum_diffusivity=diffusion
    ok=.true.
  end subroutine

  pure elemental real(dp) function les_face_mean(a,b,harmonic) result(value)
    real(dp), intent(in) :: a,b
    logical, intent(in) :: harmonic
    real(dp) :: small,large
    value=.5_dp*a+.5_dp*b
    if (.not.harmonic) return
    value=0
    if (.not.((a>0.and.b>0).or.(a<0.and.b<0))) return
    small=min(abs(a),abs(b));large=max(abs(a),abs(b))
    ! Equivalent to 2*a*b/(a+b), without product overflow/underflow.
    value=sign(small/(.5_dp+.5_dp*(small/large)),a)
  end function

  subroutine advance_dynamic_les_ssprk2_3d(species,state,temperature,spacing,test_model,test_ratio,dt,ok,grid_ratio,harmonic_mean)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp), intent(in) :: spacing(3),dt
    character(len=*), intent(in) :: test_model
    integer, intent(in) :: test_ratio
    integer, optional, intent(in) :: grid_ratio
    logical, optional, intent(in) :: harmonic_mean
    logical, intent(out) :: ok
    real(dp), allocatable :: rhs(:,:,:,:),stage(:,:,:,:),candidate(:,:,:,:),t1(:,:,:),t2(:,:,:)
    real(dp) :: diffusion
    integer :: n(3)
    logical :: local_ok
    ok=.false.;n=shape(temperature)
    if (.not.ieee_is_finite(dt)) return
    if (dt<0.or.dt>1e6_dp) return
    allocate(rhs,mold=state);allocate(stage,mold=state);allocate(candidate,mold=state)
    allocate(t1,mold=temperature);allocate(t2,mold=temperature)
    call compute_dynamic_les_rhs_3d(species,state,temperature,spacing,test_model,test_ratio,rhs,local_ok, &
      grid_ratio,diffusion, &
      harmonic_mean=harmonic_mean)
    if (.not.local_ok) return
    if (dt*diffusion*sum(1/spacing**2)>.2_dp) return
    if (dt==0) then
      ok=.true.;return
    end if
    stage=state+dt*rhs
    call recover_reactive_temperatures_3d(species,stage,temperature,n(1),n(2),n(3),t1,local_ok)
    if (.not.local_ok) return
    call compute_dynamic_les_rhs_3d(species,stage,t1,spacing,test_model,test_ratio,rhs,local_ok,grid_ratio,diffusion, &
      harmonic_mean=harmonic_mean)
    if (.not.local_ok) return
    if (dt*diffusion*sum(1/spacing**2)>.2_dp) return
    candidate=.5_dp*state+.5_dp*stage+.5_dp*dt*rhs
    call recover_reactive_temperatures_3d(species,candidate,t1,n(1),n(2),n(3),t2,local_ok)
    if (.not.local_ok) return
    state=candidate;temperature=t2;ok=.true.
  end subroutine

  ! Periodic dynamic SGS flux divergence, arithmetic cell-to-face interpolation.
  ! Pinned model transports momentum/energy only; species and mass RHS are zero.
  subroutine compute_dynamic_les_rhs_3d(species,state,temperature,spacing,test_model,test_ratio,rhs,ok, &
      grid_ratio,maximum_diffusivity,supplied_coefficients,active_x,harmonic_mean,fx,fy,fz)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),spacing(3)
    character(len=*), intent(in) :: test_model
    integer, intent(in) :: test_ratio
    real(dp), intent(inout) :: rhs(:,:,:,:)
    logical, optional, intent(in) :: harmonic_mean
    logical, intent(out) :: ok
    integer, optional, intent(in) :: grid_ratio
    real(dp), optional, intent(out) :: maximum_diffusivity
    real(dp), optional, intent(in) :: supplied_coefficients(:,:,:,:)
    integer, optional, intent(in) :: active_x(2)
    ! Full periodic face arrays, low face first, with duplicated periodic end faces.
    ! Request all three together; cropped MPI work arrays do not define full faces.
    real(dp), optional, intent(inout) :: fx(:,:,:,:),fy(:,:,:,:),fz(:,:,:,:)
    real(dp), allocatable :: faces(:,:,:,:,:)
    real(dp), allocatable :: checked_weights(:)
    integer :: first_x,last_x,face_first
    real(dp), allocatable :: coeff(:,:,:,:),q(:,:,:,:),temps(:,:,:),aij(:,:,:,:,:),alpha(:,:,:),thermal(:,:,:,:)
    real(dp), allocatable :: candidate(:,:,:,:),face_flux(:,:,:,:),viscosity(:,:,:)
    real(dp) :: diffusion,face_viscosity,delta,strain(3,3)
    real(dp) :: primitive(reactive_nprim(size(species))),g(3,3),gt(3),kij(3,3),rut(3),stress(3,3),a,ft(3)
    real(dp) :: c(3),sigma(3),velocity(3),flux(5),y(size(species)),t,sound,mw,r,cp,cv,gamma,h,e,entropy
    integer :: n(3),i,j,k,d,cell(3),lo(3),hi(3),neighbor(3),ratio
    logical :: local_ok,harmonic
    ok=.false.;n=shape(temperature);ratio=1;diffusion=0
    harmonic=.false.
    if (present(harmonic_mean)) harmonic=harmonic_mean
    if (present(maximum_diffusivity)) maximum_diffusivity=0
    if (present(grid_ratio)) ratio=grid_ratio
    if (any(shape(rhs)/=shape(state))) return
    if (present(fx).or.present(fy).or.present(fz)) then
      if (.not.(present(fx).and.present(fy).and.present(fz))) return
      if (present(active_x)) return
      if (any(shape(fx)/=[size(state,1),n(1)+1,n(2),n(3)])) return
      if (any(shape(fy)/=[size(state,1),n(1),n(2)+1,n(3)])) return
      if (any(shape(fz)/=[size(state,1),n(1),n(2),n(3)+1])) return
      allocate(faces(5,n(1),n(2),n(3),3))
    end if
    first_x=1;last_x=n(1)
    if (present(active_x)) then
      if (.not.present(supplied_coefficients)) return
      if (active_x(1)<3.or.active_x(2)>n(1)-2.or.active_x(2)<active_x(1)) return
      first_x=active_x(1);last_x=active_x(2)
    end if
    allocate(coeff(3,n(1),n(2),n(3)))
    if (present(supplied_coefficients)) then
      if (any(n<2).or.ratio<1.or.ratio>64) return
      if (size(state,1)/=reactive_nvar(size(species)).or.size(species)<1) return
      if (any([size(state,2),size(state,3),size(state,4)]/=n)) return
      if (any(shape(supplied_coefficients)/=shape(coeff))) return
      if (.not.all(ieee_is_finite(spacing)).or..not.all(ieee_is_finite(supplied_coefficients))) return
      if (any(spacing<1e-12_dp).or.any(spacing>1e9_dp)) return
      if (any(supplied_coefficients<0).or.any(supplied_coefficients>1e100_dp)) return
      if (any(supplied_coefficients(2,:,:,:)>1)) return
      call les_filter_weights(test_model,test_ratio,checked_weights,local_ok)
      if (.not.local_ok) return
      coeff=supplied_coefficients
    else
      call compute_dynamic_les_coefficients_3d(species,state,temperature,spacing,test_model,test_ratio, &
        coeff,local_ok,ratio)
      if (.not.local_ok) return
    end if
    allocate(q(reactive_nprim(size(species)),n(1),n(2),n(3)),temps(n(1),n(2),n(3)))
    allocate(aij(3,3,n(1),n(2),n(3)),alpha(n(1),n(2),n(3)),thermal(3,n(1),n(2),n(3)))
    allocate(candidate(size(state,1),n(1),n(2),n(3)),source=0._dp)
    allocate(face_flux(5,n(1),n(2),n(3)),viscosity(n(1),n(2),n(3)))
    delta=ratio*product(spacing)**(1._dp/3)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          call reactive_conserved_to_primitive(species,state(:,i,j,k),temperature(i,j,k),primitive,t,sound,local_ok)
          if (.not.local_ok) return
          if (maxval(abs(primitive(2:4)))>1e30_dp.or.primitive(1)>1e8_dp.or.t>1e8_dp) return
          q(:,i,j,k)=primitive;temps(i,j,k)=t
        end do
      end do
    end do
    do k=1,n(3)
      do j=1,n(2)
        do i=max(1,first_x-1),min(n(1),last_x+1)
          cell=[i,j,k]
          do d=1,3
            lo=cell;hi=cell;lo(d)=modulo(cell(d)-2,n(d))+1;hi(d)=modulo(cell(d),n(d))+1
            g(:,d)=(q(2:4,hi(1),hi(2),hi(3))-q(2:4,lo(1),lo(2),lo(3)))/(2*spacing(d))
            gt(d)=(temps(hi(1),hi(2),hi(3))-temps(lo(1),lo(2),lo(3)))/(2*spacing(d))
          end do
          call dynamic_smagorinsky_quantities(g,gt,q(1,i,j,k),q(2:4,i,j,k),temps(i,j,k),spacing,ratio, &
            kij,rut,stress,a,ft,local_ok)
          if (.not.local_ok) return
          aij(:,:,i,j,k)=stress;alpha(i,j,k)=a;thermal(:,i,j,k)=ft
          strain=.5_dp*(g+transpose(g))
          viscosity(i,j,k)=q(1,i,j,k)*delta**2*sqrt(2*sum(strain**2))
        end do
      end do
    end do
    do d=1,3
      face_first=first_x
      if (d==1.and.present(active_x)) face_first=first_x-1
      do k=1,n(3)
        do j=1,n(2)
          do i=face_first,last_x
            neighbor=[i,j,k];neighbor(d)=modulo(neighbor(d),n(d))+1
            c=les_face_mean(coeff(:,i,j,k),coeff(:,neighbor(1),neighbor(2),neighbor(3)),harmonic)
            stress=les_face_mean(aij(:,:,i,j,k),aij(:,:,neighbor(1),neighbor(2),neighbor(3)),harmonic)
            a=les_face_mean(alpha(i,j,k),alpha(neighbor(1),neighbor(2),neighbor(3)),harmonic)
            ft=les_face_mean(thermal(:,i,j,k),thermal(:,neighbor(1),neighbor(2),neighbor(3)),harmonic)
            velocity=.5_dp*(q(2:4,i,j,k)+q(2:4,neighbor(1),neighbor(2),neighbor(3)))
            y=.5_dp*(q(6:,i,j,k)+q(6:,neighbor(1),neighbor(2),neighbor(3)))
            t=.5_dp*(temps(i,j,k)+temps(neighbor(1),neighbor(2),neighbor(3)))
            call mixture_mass_properties(species,y,t,mw,r,cp,cv,gamma,h,e,entropy,local_ok)
            if (.not.local_ok) return
            ! Bound products before evaluation under floating-point traps.
            if (maxval(abs(c))>1e100_dp.or.maxval(abs(stress))>1e100_dp.or.abs(a)>1e100_dp) return
            if (maxval(abs(ft))>1e100_dp.or.cp>1e30_dp) return
            if (cv<1e-30_dp) return
            face_viscosity=.5_dp*(viscosity(i,j,k)+viscosity(neighbor(1),neighbor(2),neighbor(3)))
            ! Frozen-coefficient diffusion estimate including isotropic stress.
            diffusion=max(diffusion,face_viscosity*max(4*(c(1)+c(2)),cp/cv*c(3))/ &
              min(q(1,i,j,k),q(1,neighbor(1),neighbor(2),neighbor(3))))
            sigma=c(1)*stress(:,d);sigma(d)=sigma(d)-c(2)*a/3
            flux=0;flux(2:4)=-sigma;flux(5)=-dot_product(velocity,sigma)-cp*c(3)*ft(d)
            face_flux(:,i,j,k)=flux
          end do
        end do
      end do
      if (allocated(faces)) faces(:,:,:,:,d)=face_flux
      do k=1,n(3)
        do j=1,n(2)
          do i=first_x,last_x
            neighbor=[i,j,k];neighbor(d)=modulo(neighbor(d)-2,n(d))+1
            candidate(1:5,i,j,k)=candidate(1:5,i,j,k)+ &
              (face_flux(:,neighbor(1),neighbor(2),neighbor(3))-face_flux(:,i,j,k))/spacing(d)
          end do
        end do
      end do
    end do
    if (.not.all(ieee_is_finite(candidate))) return
    if (.not.ieee_is_finite(diffusion)) return
    if (allocated(faces)) then
      if (.not.all(ieee_is_finite(faces))) return
      fx=0;fy=0;fz=0
      fx(1:5,2:,:,:)=faces(:,:,:,:,1);fx(1:5,1,:,:)=faces(:,n(1),:,:,1)
      fy(1:5,:,2:,:)=faces(:,:,:,:,2);fy(1:5,:,1,:)=faces(:,:,n(2),:,2)
      fz(1:5,:,:,2:)=faces(:,:,:,:,3);fz(1:5,:,:,1)=faces(:,:,:,n(3),3)
    end if
    if (present(maximum_diffusivity)) maximum_diffusivity=diffusion
    rhs=candidate;ok=.true.
  end subroutine

  ! Periodic regular grid. Filter conserved state before EOS recovery, not T or u.
  ! Coefficients are smoothed with the pinned box filter of ratio six.
  subroutine compute_dynamic_les_coefficients_3d(species,state,temperature,spacing,test_model,test_ratio, &
      coefficients,ok,grid_ratio)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),spacing(3)
    character(len=*), intent(in) :: test_model
    integer, intent(in) :: test_ratio
    real(dp), intent(inout) :: coefficients(:,:,:,:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: grid_ratio
    real(dp), allocatable :: weights(:),smooth(:),filtered_state(:,:,:,:),q(:,:,:,:),qt(:,:,:,:)
    real(dp), allocatable :: quantities(:,:,:,:),filtered_quantities(:,:,:,:),test_quantities(:,:,:,:),raw(:,:,:,:)
    real(dp) :: primitive(reactive_nprim(size(species))),t,sound,g(3,3),gt(3),kij(3,3),aij(3,3),rut(3),a,ft(3)
    integer :: n(3),i,j,k,d,lo(3),hi(3),cell(3),ratio,level
    logical :: local_ok
    ok=.false.;n=shape(temperature);ratio=1
    if (present(grid_ratio)) ratio=grid_ratio
    if (any(n<2).or.ratio<1.or.ratio>64) return
    if (size(state,1)/=reactive_nvar(size(species)).or.size(species)<1) return
    if (any([size(state,2),size(state,3),size(state,4)]/=n)) return
    if (any(shape(coefficients)/=[3,n])) return
    if (.not.all(ieee_is_finite(spacing))) return
    if (any(spacing<1e-12_dp).or.any(spacing>1e9_dp)) return
    call les_filter_weights(test_model,test_ratio,weights,local_ok)
    if (.not.local_ok) return
    call les_filter_weights('box',6,smooth,local_ok)
    if (.not.local_ok) return
    allocate(filtered_state,mold=state)
    call filter_les_field_3d(state,weights,filtered_state,local_ok)
    if (.not.local_ok) return
    allocate(q(5,n(1),n(2),n(3)),qt(5,n(1),n(2),n(3)))
    allocate(quantities(25,n(1),n(2),n(3)),filtered_quantities(25,n(1),n(2),n(3)))
    allocate(test_quantities(25,n(1),n(2),n(3)),raw(3,n(1),n(2),n(3)))
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          call reactive_conserved_to_primitive(species,state(:,i,j,k),temperature(i,j,k),primitive,t,sound,local_ok)
          if (.not.local_ok) return
          if (maxval(abs(primitive(2:4)))>1e30_dp.or.primitive(1)>1e8_dp.or.t>1e8_dp) return
          q(:,i,j,k)=[primitive(1:4),t]
          call reactive_conserved_to_primitive(species,filtered_state(:,i,j,k),temperature(i,j,k), &
            primitive,t,sound,local_ok)
          if (.not.local_ok) return
          if (maxval(abs(primitive(2:4)))>1e30_dp.or.primitive(1)>1e8_dp.or.t>1e8_dp) return
          qt(:,i,j,k)=[primitive(1:4),t]
        end do
      end do
    end do
    do level=1,2
      if (level==2) then
        q=qt;ratio=test_ratio
      end if
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            cell=[i,j,k]
            do d=1,3
              lo=cell;hi=cell;lo(d)=modulo(cell(d)-2,n(d))+1;hi(d)=modulo(cell(d),n(d))+1
              g(:,d)=(q(2:4,hi(1),hi(2),hi(3))-q(2:4,lo(1),lo(2),lo(3)))/ (2*spacing(d))
              gt(d)=(q(5,hi(1),hi(2),hi(3))-q(5,lo(1),lo(2),lo(3)))/(2*spacing(d))
            end do
            call dynamic_smagorinsky_quantities(g,gt,q(1,i,j,k),q(2:4,i,j,k),q(5,i,j,k),spacing,ratio, &
              kij,rut,aij,a,ft,local_ok)
            if (.not.local_ok) return
            if (level==1) then
              quantities(:,i,j,k)=[reshape(kij,[9]),rut,reshape(aij,[9]),a,ft]
            else
              test_quantities(:,i,j,k)=[reshape(kij,[9]),rut,reshape(aij,[9]),a,ft]
            end if
          end do
        end do
      end do
    end do
    call filter_les_field_3d(quantities,weights,filtered_quantities,local_ok)
    if (.not.local_ok) return
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          call dynamic_smagorinsky_coefficients(qt(1,i,j,k),qt(2:4,i,j,k),qt(5,i,j,k), &
            reshape(filtered_quantities(1:9,i,j,k),[3,3]),filtered_quantities(10:12,i,j,k), &
            reshape(filtered_quantities(13:21,i,j,k),[3,3]),filtered_quantities(22,i,j,k), &
            filtered_quantities(23:25,i,j,k),reshape(test_quantities(13:21,i,j,k),[3,3]), &
            test_quantities(22,i,j,k),test_quantities(23:25,i,j,k),raw(:,i,j,k),local_ok)
          if (.not.local_ok) return
        end do
      end do
    end do
    call filter_les_field_3d(raw,smooth,coefficients,ok)
  end subroutine

  ! Supplied cell-centered derivatives: gradient(velocity component, derivative).
  ! Filter Kij, RUT, alphaij, alpha and flux_T before coefficient evaluation.
  subroutine dynamic_smagorinsky_quantities(gradient,temperature_gradient,density,velocity,temperature, &
      spacing,filter_ratio,kij,rut,alphaij,alpha,flux_t,ok)
    real(dp), intent(in) :: gradient(3,3),temperature_gradient(3),density,velocity(3),temperature,spacing(3)
    integer, intent(in) :: filter_ratio
    real(dp), intent(out) :: kij(3,3),rut(3),alphaij(3,3),alpha,flux_t(3)
    logical, intent(out) :: ok
    real(dp) :: viscosity,strain(3,3),stress(3,3),magnitude
    integer :: i,j
    kij=0;rut=0;alphaij=0;alpha=0;flux_t=0;ok=.false.
    if (.not.all(ieee_is_finite([temperature_gradient,velocity,temperature]))) return
    if (maxval(abs(temperature_gradient))>1e30_dp.or.maxval(abs(velocity))>1e30_dp) return
    if (temperature<=0.or.temperature>1e8_dp.or.filter_ratio<1.or.filter_ratio>64) return
    call smagorinsky_stress(gradient,density,spacing,1._dp,viscosity,stress,ok)
    if (.not.ok) return
    viscosity=viscosity*real(filter_ratio,dp)**2
    alphaij=stress*real(filter_ratio,dp)**2
    strain=.5_dp*(gradient+transpose(gradient));magnitude=sqrt(2*sum(strain**2))
    alpha=2*viscosity*magnitude;flux_t=viscosity*temperature_gradient
    do j=1,3
      do i=1,3
        kij(i,j)=density*velocity(i)*velocity(j)
      end do
    end do
    rut=density*velocity*temperature
  end subroutine

  ! Inputs alpha* and flux_T are test-filtered grid quantities; beta* and
  ! test_flux_T come from the test-filtered primitive state at the test width.
  ! Output ordering: Cs^2, CI, Cs^2/PrT. This routine performs no spatial filter.
  subroutine dynamic_smagorinsky_coefficients(density,velocity,temperature,kij,rut,alphaij,alpha,flux_t, &
      betaij,beta,test_flux_t,coefficients,ok)
    real(dp), intent(in) :: density,velocity(3),temperature,kij(3,3),rut(3),alphaij(3,3),alpha,flux_t(3)
    real(dp), intent(in) :: betaij(3,3),beta,test_flux_t(3)
    real(dp), intent(out) :: coefficients(3)
    logical, intent(out) :: ok
    real(dp) :: l(3,3),m(3,3),t(3),ke(3),lkk,bma,candidate(3)
    integer :: i,j
    coefficients=0;ok=.false.
    if (.not.all(ieee_is_finite([density,velocity,temperature,rut,alpha,flux_t,beta,test_flux_t]))) return
    if (.not.all(ieee_is_finite(kij)).or..not.all(ieee_is_finite(alphaij))) return
    if (.not.all(ieee_is_finite(betaij))) return
    if (density<=0.or.density>1e8_dp.or.temperature<=0.or.temperature>1e8_dp) return
    if (maxval(abs(velocity))>1e30_dp) return
    ! Bound contractions before evaluating products under floating-point traps.
    if (maxval(abs([rut,alpha,flux_t,beta,test_flux_t]))>1e100_dp) return
    if (maxval(abs(kij))>1e100_dp.or.maxval(abs(alphaij))>1e100_dp.or.maxval(abs(betaij))>1e100_dp) return
    do j=1,3
      do i=1,3
        l(i,j)=kij(i,j)-density*velocity(i)*velocity(j)
      end do
    end do
    m=betaij-alphaij;t=flux_t-test_flux_t;ke=rut-density*velocity*temperature
    lkk=l(1,1)+l(2,2)+l(3,3);bma=beta-alpha+regularizer
    if (bma==0) return
    candidate(1)=max(sum(l*m)/(sum(m*m)+regularizer),regularizer)
    ! Apply the pinned [floor,1] clipping before division to avoid overflow
    ! when the regularized isotropic denominator is nonzero but very small.
    candidate(2)=regularizer
    if ((lkk>0.and.bma>0).or.(lkk<0.and.bma<0)) then
      candidate(2)=1
      if (abs(lkk)<abs(bma)) candidate(2)=max(lkk/bma,regularizer)
    end if
    candidate(3)=max(sum(ke*t)/(sum(t*t)+regularizer),regularizer)
    if (.not.all(ieee_is_finite(candidate))) return
    coefficients=candidate;ok=.true.
  end subroutine
end module
