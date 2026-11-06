! Independent implementation of pinned PeleC pc_rst_int_e policy.
module dual_energy_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use mixture_thermo_mod, only: mixture_mass_properties
  use nasa7_thermo_mod, only: nasa7_species,valid_nasa7_species,universal_gas_constant
  implicit none
  private
  type, public :: dual_energy_policy
    logical :: allow_small_energy=.true.
    logical :: allow_negative_energy=.true.
    logical :: update_total_from_internal=.true.
    real(dp) :: eta2=1.e-4_dp
  end type
  abstract interface
    subroutine dual_energy_evaluator(temperature,mass_fractions,energy,ok)
      import :: dp
      real(dp), intent(in) :: temperature,mass_fractions(:)
      real(dp), intent(out) :: energy
      logical, intent(out) :: ok
    end subroutine
  end interface
  public :: dual_energy_evaluator,reset_dual_energy_cell,reset_dual_energy_cells,reset_nasa7_dual_energy_cell, &
    reset_nasa7_dual_energy_patch,initialize_nasa7_internal_energy_patch
contains
  ! Native UEINT initialization uses rho*e(T,Y), independently of kinetic
  ! and rotating total energy. Covered values and all outputs survive failure.
  subroutine initialize_nasa7_internal_energy_patch(species,state,temperature,fluid_volume,rho_internal,ok)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),fluid_volume(:,:,:)
    real(dp), intent(inout) :: rho_internal(:,:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: candidate(:,:,:)
    real(dp) :: rho,y(size(species)),mw,r,cp,cv,gamma,h,e,entropy
    integer :: n(3),i,j,k
    logical :: valid
    ok=.false.;n=shape(rho_internal)
    if(size(species)<1.or.any(n<1)) return
    if(any(shape(state)/=[5+size(species),n])) return
    if(any(shape(temperature)/=n).or.any(shape(fluid_volume)/=n)) return
    if(any(.not.ieee_is_finite(fluid_volume)).or.any(fluid_volume<0)) return
    candidate=rho_internal
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      if(fluid_volume(i,j,k)==0) cycle
      if(any(.not.ieee_is_finite(state(:,i,j,k)))) return
      rho=state(1,i,j,k)
      if(rho<=0.or.any(state(6:,i,j,k)<0)) return
      ! Bound ratios before division and reduction, including overflow-trap
      ! builds. Invalid finite inputs must refuse before publishing a prefix.
      if(rho<1) then
        if(any(state(6:,i,j,k)>huge(1._dp)*rho)) return
      endif
      y=state(6:,i,j,k)/rho
      if(any(y>1+5.e-12_dp)) return
      if(abs(sum(y)-1)>5.e-12_dp) return
      call mixture_mass_properties(species,y,temperature(i,j,k),mw,r,cp,cv,gamma,h,e,entropy,valid)
      if(.not.valid) return
      if(rho>huge(1._dp)/max(1._dp,abs(e))) return
      candidate(i,j,k)=rho*e
      if(.not.ieee_is_finite(candidate(i,j,k))) return
    enddo
    enddo
    enddo
    rho_internal=candidate;ok=.true.
  end subroutine initialize_nasa7_internal_energy_patch

  subroutine reset_nasa7_dual_energy_cell(species,state,rho_internal,temperature,policy,energy_change,ok,potential)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:),rho_internal,energy_change
    real(dp), intent(in) :: temperature
    type(dual_energy_policy), intent(in) :: policy
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: potential
    integer :: sp
    ok=.false.
    if(size(state)/=5+size(species)) return
    do sp=1,size(species)
      if(.not.valid_nasa7_species(species(sp))) return
    enddo
    call reset_dual_energy_cell(state,rho_internal,temperature,policy,energy_change,ok,nasa_energy,potential)
  contains
    subroutine nasa_energy(t,y,e,success)
      real(dp), intent(in) :: t,y(:)
      real(dp), intent(out) :: e
      logical, intent(out) :: success
      real(dp) :: a(7),polynomial,term,next,gas_constant,limit
      integer :: sp,power
      success=.false.;e=0;limit=huge(1._dp)/4
      ! T2Ei uses NASA7 coefficients outside their fitted temperature range
      ! for this explicit reset policy. Evaluate e(T) directly: h/T and
      ! entropy are unnecessary and ill-conditioned near machine tiny.
      do sp=1,size(species)
        a=species(sp)%low_coefficients
        if(t>=species(sp)%temperature_mid) a=species(sp)%high_coefficients
        polynomial=a(5)/5
        do power=4,1,-1
          if(abs(polynomial)>limit/max(1._dp,t)) return
          term=a(power)/real(power,dp)
          if(power==1) term=term-1
          next=polynomial*t
          if(abs(next)>limit-abs(term)) return
          polynomial=next+term
        enddo
        if(abs(polynomial)>limit/max(1._dp,t)) return
        next=polynomial*t
        if(abs(next)>limit-abs(a(6))) return
        polynomial=next+a(6)
        gas_constant=universal_gas_constant/species(sp)%molecular_weight
        if(abs(polynomial)>limit/max(1._dp,abs(gas_constant))) return
        term=y(sp)*gas_constant*polynomial
        if(abs(e)>limit-abs(term)) return
        e=e+term
      enddo
      success=ieee_is_finite(e)
    end subroutine nasa_energy
  end subroutine reset_nasa7_dual_energy_cell

  ! Reset only fluid cells. Cell volumes include the EB volume fraction; all
  ! outputs commit together so covered-cell sentinels survive a failed stage.
  subroutine reset_nasa7_dual_energy_patch(species,state,rho_internal,temperature,fluid_volume,policy, &
      energy_change,ok,potential)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(inout) :: state(:,:,:,:),rho_internal(:,:,:),energy_change
    real(dp), intent(in) :: temperature(:,:,:),fluid_volume(:,:,:)
    type(dual_energy_policy), intent(in) :: policy
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: potential(:,:,:)
    real(dp), allocatable :: candidate(:,:,:,:),held(:,:,:)
    real(dp) :: correction,delta,phi
    integer :: n(3),i,j,k
    logical :: valid
    ok=.false.;n=shape(rho_internal)
    if(any(shape(state)/=[5+size(species),n])) return
    if(any(shape(temperature)/=n).or.any(shape(fluid_volume)/=n)) return
    if(present(potential)) then
      if(any(shape(potential)/=n)) return
    endif
    if(.not.ieee_is_finite(energy_change)) return
    if(any(.not.ieee_is_finite(fluid_volume)).or.any(fluid_volume<0)) return
    candidate=state;held=rho_internal;correction=0
    do k=1,n(3)
    do j=1,n(2)
    do i=1,n(1)
      if(fluid_volume(i,j,k)==0) cycle
      phi=0
      if(present(potential)) phi=potential(i,j,k)
      delta=0
      call reset_nasa7_dual_energy_cell(species,candidate(:,i,j,k),held(i,j,k),temperature(i,j,k), &
        policy,delta,valid,phi)
      if(.not.valid) return
      correction=correction+fluid_volume(i,j,k)*delta
      if(.not.ieee_is_finite(correction)) return
    enddo
    enddo
    enddo
    if(.not.ieee_is_finite(energy_change+correction)) return
    state=candidate;rho_internal=held;energy_change=energy_change+correction;ok=.true.
  end subroutine reset_nasa7_dual_energy_patch

  subroutine reset_dual_energy_cell(state,rho_internal,temperature,policy,energy_change,ok, &
      evaluate_energy,potential)
    real(dp), intent(inout) :: state(:),rho_internal,energy_change
    real(dp), intent(in) :: temperature
    type(dual_energy_policy), intent(in) :: policy
    logical, intent(out) :: ok
    procedure(dual_energy_evaluator), optional :: evaluate_energy
    real(dp), optional, intent(in) :: potential
    real(dp) :: rho,velocity(3),ke,phi,total,internal,from_total,rho_floor,small_e,replacement,change
    real(dp) :: fractions(max(0,size(state)-5))
    logical :: valid
    ok=.false.
    if(size(state)<6) return
    if(.not.all(ieee_is_finite(state)).or..not.ieee_is_finite(rho_internal)) return
    if(.not.all(ieee_is_finite([temperature,policy%eta2]))) return
    rho=state(1)
    if(rho<=0.or.any(state(6:)<0)) return
    if(abs(sum(state(6:))-rho)>5.e-12_dp*rho) return
    phi=0
    if(present(potential)) phi=potential
    if(.not.ieee_is_finite(phi)) return
    velocity=state(2:4)/rho
    ke=.5_dp*sum(velocity**2)-phi
    total=state(5);internal=rho_internal;from_total=total-rho*ke
    if(.not.all(ieee_is_finite([ke,from_total]))) return
    rho_floor=rho*tiny(1._dp)
    fractions=state(6:)/rho
    if(.not.policy%allow_small_energy) then
      call energy_at(epsilon(1._dp)*1.e-100_dp,small_e,valid)
      if(.not.valid) return
      if(total/rho<small_e) then
        if(internal/rho<small_e) then
          call energy_at(max(temperature,tiny(1._dp)),replacement,valid)
          if(.not.valid) return
          internal=rho*replacement
        endif
        ! Pinned branch repairs total energy even when update flag is false.
        total=internal+rho*ke
      else
        if(trust_total(from_total,total,rho_floor,policy%eta2)) internal=from_total
        if(internal/rho<small_e) then
          call energy_at(max(temperature,tiny(1._dp)),replacement,valid)
          if(.not.valid) return
          if(policy%update_total_from_internal) total=total+(rho*replacement-internal)
          internal=rho*replacement
        endif
      endif
    else if(.not.policy%allow_negative_energy) then
      if(total<rho_floor) then
        if(internal<rho_floor) then
          call energy_at(epsilon(1._dp)*1.e-100_dp,replacement,valid)
          if(.not.valid) return
          internal=rho*replacement
        endif
        total=internal+rho*ke
      else
        if(trust_total(from_total,total,rho_floor,policy%eta2)) then
          internal=from_total
        else if(internal>rho_floor.and.policy%update_total_from_internal) then
          total=internal+rho*ke
        else if(internal<=rho_floor) then
          call energy_at(epsilon(1._dp)*1.e-100_dp,replacement,valid)
          if(.not.valid) return
          if(policy%update_total_from_internal) total=total+(rho*replacement-internal)
          internal=rho*replacement
        endif
      endif
    else
      internal=from_total
    endif
    change=total-state(5)
    if(.not.all(ieee_is_finite([total,internal,change]))) return
    state(5)=total;rho_internal=internal;energy_change=change;ok=.true.
  contains
    subroutine energy_at(t,e,success)
      real(dp), intent(in) :: t
      real(dp), intent(out) :: e
      logical, intent(out) :: success
      e=0;success=.false.
      if(.not.present(evaluate_energy)) return
      call evaluate_energy(t,fractions,e,success)
      success=success.and.ieee_is_finite(e)
    end subroutine
  end subroutine reset_dual_energy_cell

  logical function trust_total(internal,total,floor,eta) result(trusted)
    real(dp), intent(in) :: internal,total,floor,eta
    trusted=.false.
    if(internal<=floor) return
    ! The pinned quotient is +infinity at a zero denominator with positive
    ! numerator. Select its branch without raising a divide-by-zero trap.
    if(total==floor) then
      trusted=.true.
    else
      trusted=(internal-floor)/(total-floor)>eta
    endif
  end function trust_total

  subroutine reset_dual_energy_cells(state,rho_internal,temperature,policy,energy_change,ok,evaluate_energy,potential)
    real(dp), intent(inout) :: state(:,:),rho_internal(:),energy_change(:)
    real(dp), intent(in) :: temperature(:)
    type(dual_energy_policy), intent(in) :: policy
    logical, intent(out) :: ok
    procedure(dual_energy_evaluator), optional :: evaluate_energy
    real(dp), optional, intent(in) :: potential(:)
    real(dp), allocatable :: candidate(:,:),aux(:),change(:)
    real(dp) :: phi
    integer :: cell,count
    logical :: valid
    ok=.false.;count=size(state,2)
    if(size(rho_internal)/=count.or.size(temperature)/=count.or.size(energy_change)/=count) return
    if(present(potential)) then
      if(size(potential)/=count) return
    endif
    candidate=state;aux=rho_internal;allocate(change(count))
    do cell=1,count
      phi=0
      if(present(potential)) phi=potential(cell)
      call reset_dual_energy_cell(candidate(:,cell),aux(cell),temperature(cell),policy,change(cell),valid, &
        evaluate_energy,phi)
      if(.not.valid) return
    enddo
    state=candidate;rho_internal=aux;energy_change=change;ok=.true.
  end subroutine reset_dual_energy_cells
end module dual_energy_mod
