module spray_reactive_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use reactive_3d_mod, only: advance_reactive_chemistry_3d, advance_reactive_full_3d, recover_reactive_temperatures_3d
  use spray_parcel_mod
  use spray_coupling_3d_mod, only: advance_spray_coupling_3d
  use les_smagorinsky_3d_mod, only: advance_les_ssprk2_3d
#ifdef PELEF_SPRAY_CVODE
  use sundials_constant_volume_reactor_mod, only: cvode_reactor_context, initialize_constant_volume_cvode, &
    advance_constant_volume_cvode, finalize_constant_volume_cvode
#endif
  implicit none
  private
  public :: advance_spray_chemistry_3d, advance_spray_reactive_3d, spray_cvode_available
contains
  logical function spray_cvode_available() result(available)
#ifdef PELEF_SPRAY_CVODE
    available=.true.
#else
    available=.false.
#endif
  end function

  subroutine advance_spray_chemistry_3d(species,reactions,state,temperature,interval,rtol,atol,method,ok,message)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    real(dp), intent(in) :: interval,rtol,atol
    character(len=*), intent(in) :: method
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    integer :: dims(3)
#ifdef PELEF_SPRAY_CVODE
    type(cvode_reactor_context) :: context
    real(dp), allocatable :: candidate(:,:,:,:),temp(:,:,:),check_t(:,:,:)
    real(dp) :: primitive(reactive_nprim(size(species))),y(size(species)),rho,e,t,sound,initial_step
    integer :: i,j,k
    logical :: local_ok,finalize_ok
    character(len=256) :: finalize_message
#endif
    ok=.false.; message='Invalid spray chemistry request'
    dims=shape(temperature)
    if (size(state,1) /= reactive_nvar(size(species))) return
    if (any(shape(state(1,:,:,:)) /= dims) .or. any(dims < 1)) return
    if (.not. all(ieee_is_finite([interval,rtol,atol]))) return
    if (interval < 0.0_dp .or. interval > 1.0e6_dp .or. min(rtol,atol) <= 0.0_dp) return
    select case(trim(method))
    case('implicit','explicit')
      call advance_reactive_chemistry_3d(species,reactions,state,temperature,dims(1),dims(2),dims(3), &
        interval,rtol,atol,ok,chemistry_integrator=method)
      if (ok) message=''
    case('cvode')
#ifndef PELEF_SPRAY_CVODE
      message='CVODE requested but this spray runtime was built without SUNDIALS'
#else
      allocate(candidate,source=state); allocate(temp,source=temperature); allocate(check_t,mold=temperature)
      call recover_reactive_temperatures_3d(species,state,temperature,dims(1),dims(2),dims(3),check_t,local_ok)
      if (.not. local_ok) return
      if (interval == 0.0_dp) then
        ok=.true.; message=''; return
      end if
      initial_step=min(interval,1.0e-12_dp)
      do k=1,dims(3)
        do j=1,dims(2)
          do i=1,dims(1)
            call reactive_conserved_to_primitive(species,state(:,i,j,k),temperature(i,j,k),primitive,t,sound,local_ok)
            if (.not. local_ok) return
            rho=primitive(1); y=primitive(6:)
            e=state(5,i,j,k)/rho-0.5_dp*sum(primitive(2:4)**2)
            ! A new context after transport/spray uses the current rho and energy.
            ! Finalize every cell, including failed advances: no stale state or slot leak.
            call initialize_constant_volume_cvode(context,species,reactions,rho,e,0.0_dp,initial_step, &
              min(initial_step,interval*1.0e-12_dp),interval,rtol,atol,100000,y,t,local_ok,message)
            if (.not. local_ok) return
            call advance_constant_volume_cvode(context,interval,y,t,local_ok,message)
            call finalize_constant_volume_cvode(context,finalize_ok,finalize_message)
            if (.not. finalize_ok) then
              message=finalize_message; return
            end if
            if (.not. local_ok) return
            candidate(6:,i,j,k)=rho*y; temp(i,j,k)=t
          end do
        end do
      end do
      call recover_reactive_temperatures_3d(species,candidate,temp,dims(1),dims(2),dims(3),check_t,local_ok)
      if (.not. local_ok) return
      ! Chemistry never overwrites density, momentum or total energy.
      state=candidate; temperature=temp; ok=.true.; message=''
#endif
    case default
      message='Unknown chemistry backend (choose implicit, explicit, or cvode)'
    end select
  end subroutine

  subroutine advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,parcels, &
      lower,spacing,dt,chemistry_enabled,transport_enabled,chemistry_method,rtol,atol,cs,prt,sct,ok,message, &
      reconstruction,riemann_solver)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    type(spray_liquid), intent(in) :: liquid
    type(spray_options), intent(in) :: options
    real(dp), intent(inout) :: state(:,:,:,:),temperature(:,:,:)
    type(spray_parcel), intent(inout) :: parcels(:)
    real(dp), intent(in) :: lower(3),spacing(3),dt,rtol,atol,cs,prt,sct
    logical, intent(in) :: chemistry_enabled,transport_enabled
    character(len=*), intent(in) :: chemistry_method
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    real(dp), allocatable :: candidate(:,:,:,:),temp(:,:,:)
    type(spray_parcel), allocatable :: drops(:)
    character(len=*), intent(in), optional :: reconstruction,riemann_solver
    character(len=32) :: selected_reconstruction,selected_riemann
    real(dp) :: theta
    integer :: dims(3),substeps,half
    logical :: local_ok
    ok=.false.; message='Invalid coupled spray step'
    selected_reconstruction='characteristic_plm'; selected_riemann='hllc'
    if (present(reconstruction)) selected_reconstruction=reconstruction
    if (present(riemann_solver)) selected_riemann=riemann_solver
    if (.not. all(ieee_is_finite([dt,rtol,atol,cs,prt,sct]))) return
    if (dt <= 0.0_dp .or. dt > 1.0e6_dp) return
    dims=shape(temperature)
    allocate(candidate,source=state); allocate(temp,source=temperature); allocate(drops,source=parcels)
    ! Symmetric operator ordering; the nonlinear frozen-film parcel integrator
    ! remains first order in its substep. No second-order spray claim is implied.
    do half=1,2
      if (half == 1) then
        call advance_spray_coupling_3d(species,transport,liquid,options,candidate,temp,lower,spacing, &
          drops,0.5_dp*dt,substeps,local_ok,message)
        if (.not. local_ok) return
        call advance_les_ssprk2_3d(species,candidate,temp,spacing,cs,prt,sct,0.5_dp*dt,local_ok)
        if (.not. local_ok) then
          message='LES step rejected'; return
        end if
      end if
      if (chemistry_enabled) then
        call advance_spray_chemistry_3d(species,reactions,candidate,temp,0.5_dp*dt,rtol,atol,chemistry_method,local_ok,message)
        if (.not. local_ok) return
      end if
      if (half == 1) then
        call advance_reactive_full_3d(species,reactions,transport,candidate,temp,dims(1),dims(2),dims(3), &
          spacing(1),spacing(2),spacing(3),dt,trim(selected_riemann),.false.,rtol,atol,transport_enabled, &
          .true.,.true.,.true.,.false.,theta,local_ok,trim(selected_reconstruction),'mc',chemistry_integrator='implicit')
        if (.not. local_ok) then
          message='Hydro/molecular transport step rejected'; return
        end if
      else
        call advance_les_ssprk2_3d(species,candidate,temp,spacing,cs,prt,sct,0.5_dp*dt,local_ok)
        if (.not. local_ok) then
          message='LES step rejected'; return
        end if
        call advance_spray_coupling_3d(species,transport,liquid,options,candidate,temp,lower,spacing, &
          drops,0.5_dp*dt,substeps,local_ok,message)
        if (.not. local_ok) return
      end if
    end do
    state=candidate; temperature=temp; parcels=drops; ok=.true.; message=''
  end subroutine
end module spray_reactive_3d_mod
