module sundials_forced_reactor_mod
  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species,nasa7_mass_properties
  use elementary_kinetics_mod, only: elementary_reaction, &
    elementary_production_rates_from_concentrations,elementary_production_jacobian
  use constant_volume_reactor_mod, only: forced_reactor_rhs,forced_reactor_reduced_jacobian
  use fsundials_core_mod
  use fcvode_mod
  use fnvector_serial_mod, only: FN_VMake_Serial
  use fsunmatrix_dense_mod, only: FSUNDenseMatrix,FSUNDenseMatrix_Data
  use fsunlinsol_dense_mod, only: FSUNLinSol_Dense
  implicit none
  private
#ifdef SUNDIALS_INT32_T
  integer, parameter :: index_kind=c_int32_t
#else
  integer, parameter :: index_kind=c_int64_t
#endif
  type :: forcing_context
    logical :: busy=.false.
    type(nasa7_species), allocatable :: species(:)
    type(elementary_reaction), allocatable :: reactions(:)
    real(dp), allocatable :: source(:)
    integer, allocatable :: independent(:)
    integer :: dependent=0
    real(dp) :: density=0,density_source=0
    real(dp) :: energy=0,energy_source=0,guess=0
  end type
  type(forcing_context), save :: registry(64)
  public :: advance_forced_reactor_cvode,advance_forced_reactor_cvode_full
  public :: evaluate_forced_reactor_full_rhs
contains
  ! One complete interval is transactional. Context slots are private to the
  ! synchronous call; callback tokens never escape the lifetime of CVODE memory.
  ! As with the existing CVODE registry, concurrent threaded calls are unsupported.
  subroutine advance_forced_reactor_cvode(species,reactions,initial_energy_density,external_species, &
      external_energy,duration,rtol,atol,species_density,temperature,completed_steps,ok,maximum_steps)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: initial_energy_density,external_species(:),external_energy,duration,rtol,atol
    real(dp), intent(inout) :: species_density(:),temperature
    integer, intent(out) :: completed_steps
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_steps
    real(c_double), allocatable, target :: state(:),constraints(:)
    real(dp) :: rates(size(species)),recovered,full_state(size(species)),final_density,final_temperature
    real(c_double) :: returned(1)
    integer(c_long) :: steps(1)
    integer(c_int), target :: slot
    integer(c_int) :: status,cleanup
    integer :: i,n,limit
    logical :: valid,accepted
    type(c_ptr) :: sunctx,memory
    type(N_Vector), pointer :: vector,limits
    type(SUNMatrix), pointer :: matrix
    type(SUNLinearSolver), pointer :: solver
    ok=.false.;completed_steps=0;limit=100000
    if(present(maximum_steps)) limit=maximum_steps
    if(.not.all(ieee_is_finite([duration,rtol,atol]))) return
    if(duration<0.or.min(rtol,atol)<=0.or.limit<1) return
    call forced_reactor_rhs(species,reactions,initial_energy_density,0._dp,species_density,external_species, &
      external_energy,temperature,rates,recovered,valid)
    if(.not.valid) return
    if(duration==0) then
      ok=.true.;return
    endif
    final_density=sum(species_density)+duration*sum(external_species)
    if(.not.ieee_is_finite(final_density).or.final_density<=0) return
    if(size(species)==1) then
      full_state=final_density
      call forced_reactor_rhs(species,reactions,initial_energy_density,duration,full_state,external_species, &
        external_energy,recovered,rates,final_temperature,valid)
      if(.not.valid) return
      species_density=full_state;temperature=final_temperature;completed_steps=1;ok=.true.;return
    end if
    slot=0
    do i=1,size(registry)
      if(registry(i)%busy) cycle
      slot=int(i,c_int);exit
    enddo
    if(slot==0) return
    registry(slot)%busy=.true.;registry(slot)%species=species;registry(slot)%reactions=reactions
    registry(slot)%source=external_species;registry(slot)%energy=initial_energy_density
    registry(slot)%energy_source=external_energy;registry(slot)%guess=recovered
    registry(slot)%density=sum(species_density);registry(slot)%density_source=sum(external_species)
    registry(slot)%dependent=maxloc(species_density,dim=1)
    registry(slot)%independent=pack([(i,i=1,size(species))],[(i/=registry(slot)%dependent,i=1,size(species))])
    n=size(species)-1;allocate(state(n),constraints(n))
    state=species_density(registry(slot)%independent);constraints=1
    sunctx=c_null_ptr;memory=c_null_ptr;nullify(vector,limits,matrix,solver);accepted=.false.
    solve: block
      status=FSUNContext_Create(SUN_COMM_NULL,sunctx)
      if(status/=0.or..not.c_associated(sunctx)) exit solve
      vector=>FN_VMake_Serial(int(n,index_kind),state,sunctx)
      limits=>FN_VMake_Serial(int(n,index_kind),constraints,sunctx)
      if(.not.associated(vector).or..not.associated(limits)) exit solve
      matrix=>FSUNDenseMatrix(int(n,index_kind),int(n,index_kind),sunctx)
      if(.not.associated(matrix)) exit solve
      solver=>FSUNLinSol_Dense(vector,matrix,sunctx)
      if(.not.associated(solver)) exit solve
      memory=FCVodeCreate(CV_BDF,sunctx)
      if(.not.c_associated(memory)) exit solve
      status=FCVodeInit(memory,c_funloc(forced_cvode_rhs),0._c_double,vector)
      if(status/=CV_SUCCESS) exit solve
      status=FCVodeSetUserData(memory,c_loc(slot))
      if(status/=CV_SUCCESS) exit solve
      status=FCVodeSStolerances(memory,real(rtol,c_double),real(atol,c_double))
      if(status/=CV_SUCCESS) exit solve
      ! A frozen SDC/MOL source can remove an initially absent species before
      ! chemistry produces it later in this interval. Such auxiliary forced
      ! trajectories need signed internal steps (as in native PelePhysics).
      ! Constraining every CVODE substep can stall at t=0 even when the final
      ! state is admissible. Keep constraints for nonnegative forcing; every
      ! endpoint still passes the same nonnegativity/EOS/rollback checks below.
      if(all(external_species>=0._dp)) then
        status=FCVodeSetConstraints(memory,limits)
        if(status/=CV_SUCCESS) exit solve
      endif
      status=FCVodeSetLinearSolver(memory,solver,matrix)
      if(status/=CV_SUCCESS) exit solve
      status=FCVodeSetJacFn(memory,c_funloc(forced_cvode_jacobian))
      if(status/=CV_SUCCESS) exit solve
      ! Evolve independent densities; recover the largest initial species from
      ! the analytic total-density source balance, as in the fixed-density backend.
      status=FCVodeSetMaxNumSteps(memory,int(limit,c_long))
      if(status/=CV_SUCCESS) exit solve
      status=FCVodeSetMaxStep(memory,real(duration,c_double))
      if(status/=CV_SUCCESS) exit solve
      status=FCVodeSetStopTime(memory,real(duration,c_double))
      if(status/=CV_SUCCESS) exit solve
      status=FCVode(memory,real(duration,c_double),vector,returned,CV_NORMAL)
      if(status/=CV_SUCCESS.and.status/=CV_TSTOP_RETURN) exit solve
      if(abs(real(returned(1),dp)-duration)>100*epsilon(1._dp)*duration) exit solve
      status=FCVodeGetNumSteps(memory,steps)
      if(status/=CV_SUCCESS.or.steps(1)>int(limit,c_long)) exit solve
      ! Constraint-satisfying solves can return signed roundoff for absent
      ! trace species. Repair only machine-scale endpoint noise, independently
      ! of user tolerances; reconstruction retains the analytic total mass.
      if(any(state < -16*epsilon(1._dp)*final_density)) exit solve
      state=max(state,0._c_double)
      call reconstruct_density(registry(slot),duration,state,full_state)
      call forced_reactor_rhs(species,reactions,initial_energy_density,duration,full_state,external_species, &
        external_energy,registry(slot)%guess,rates,recovered,valid)
      if(.not.valid) exit solve
      accepted=.true.
    end block solve
    if(c_associated(memory)) call FCVodeFree(memory)
    if(associated(solver)) then
      cleanup=FSUNLinSolFree(solver)
      if(cleanup/=0) accepted=.false.
    endif
    if(associated(matrix)) call FSUNMatDestroy(matrix)
    if(associated(limits)) call FN_VDestroy(limits)
    if(associated(vector)) call FN_VDestroy(vector)
    if(c_associated(sunctx)) then
      cleanup=FSUNContext_Free(sunctx)
      if(cleanup/=0) accepted=.false.
    endif
    deallocate(registry(slot)%species,registry(slot)%reactions,registry(slot)%source,registry(slot)%independent)
    registry(slot)%busy=.false.
    if(.not.accepted) return
    species_density=full_state;temperature=recovered;completed_steps=int(steps(1));ok=.true.
  end subroutine

  integer(c_int) function forced_cvode_rhs(time,y,dy,user_data) result(status) bind(C)
    real(c_double), value :: time
    type(N_Vector) :: y,dy
    type(c_ptr), value :: user_data
    integer(c_int), pointer :: slot
    real(c_double), pointer :: state(:),rhs(:)
    real(dp) :: recovered
    real(dp), allocatable :: full_state(:),rates(:)
    logical :: valid
    status=-1
    if(.not.c_associated(user_data)) return
    call c_f_pointer(user_data,slot)
    if(slot<1.or.slot>size(registry)) return
    if(.not.registry(slot)%busy) return
    state=>FN_VGetArrayPointer(y);rhs=>FN_VGetArrayPointer(dy)
    if(.not.associated(state).or..not.associated(rhs)) return
    allocate(full_state(size(registry(slot)%species)),rates(size(registry(slot)%species)))
    if(.not.all(ieee_is_finite(state))) return
    ! Newton trial coordinates may be negative before CVODE enforces its
    ! solution constraints. Evaluate their signed polynomial/thermodynamic state;
    ! the accepted endpoint still passes strict density/thermodynamic checks.
    call reconstruct_density(registry(slot),real(time,dp),state,full_state)
    call forced_reactor_rhs(registry(slot)%species,registry(slot)%reactions,registry(slot)%energy, &
      real(time,dp),full_state,registry(slot)%source,registry(slot)%energy_source, &
      registry(slot)%guess,rates,recovered,valid,permit_signed_trials=.true.)
    status=1
    if(.not.valid) return
    rhs=rates(registry(slot)%independent)
    registry(slot)%guess=recovered;status=0
  end function
  integer(c_int) function forced_cvode_jacobian(time,y,dy,matrix,user_data,tmp1,tmp2,tmp3) &
      result(status) bind(C)
    real(c_double), value :: time
    type(N_Vector) :: y,dy,tmp1,tmp2,tmp3
    type(SUNMatrix) :: matrix
    type(c_ptr), value :: user_data
    integer(c_int), pointer :: slot
    real(c_double), pointer :: state(:),matrix_values(:,:)
    real(dp), allocatable :: full_state(:),jacobian(:,:)
    real(dp) :: recovered
    integer :: n
    logical :: valid
    status=-1
    if(.not.c_associated(user_data)) return
    call c_f_pointer(user_data,slot)
    if(slot<1.or.slot>size(registry)) return
    if(.not.registry(slot)%busy.or..not.ieee_is_finite(real(time,dp))) return
    state=>FN_VGetArrayPointer(y)
    if(.not.associated(state)) return
    n=size(registry(slot)%independent)
    if(size(state)/=n.or..not.all(ieee_is_finite(state))) return
    allocate(full_state(n+1),jacobian(n,n))
    call reconstruct_density(registry(slot),real(time,dp),state,full_state)
    call forced_reactor_reduced_jacobian(registry(slot)%species,registry(slot)%reactions, &
      registry(slot)%energy,real(time,dp),full_state,registry(slot)%source,registry(slot)%energy_source, &
      registry(slot)%guess,jacobian,recovered,valid,registry(slot)%dependent,permit_signed_trials=.true.)
    status=1
    if(.not.valid) return
    matrix_values(1:n,1:n)=>FSUNDenseMatrix_Data(matrix)
    status=-1
    if(.not.associated(matrix_values)) return
    matrix_values=real(jacobian,c_double)
    ! The signed trial derivative uses the same unprojected coordinates as
    ! the RHS. Accepted endpoint constraints and mass admission are unchanged.
    registry(slot)%guess=recovered;status=0
    if(associated(FN_VGetArrayPointer(dy))) continue
    if(associated(FN_VGetArrayPointer(tmp1))) continue
    if(associated(FN_VGetArrayPointer(tmp2))) continue
    if(associated(FN_VGetArrayPointer(tmp3))) continue
  end function forced_cvode_jacobian
  subroutine reconstruct_density(context,time,reduced,full_state)
    type(forcing_context), intent(in) :: context
    real(dp), intent(in) :: time
    real(c_double), intent(in) :: reduced(:)
    real(dp), intent(out) :: full_state(:)
    full_state(context%independent)=real(reduced,dp)
    full_state(context%dependent)=context%density+time*context%density_source-sum(real(reduced,dp))
  end subroutine
  ! Explicit full-coordinate entry: atol is a density-space tolerance supplied by the caller.
  subroutine advance_forced_reactor_cvode_full(species,reactions,initial_energy_density,external_species, &
      external_energy,duration,rtol,atol,species_density,temperature,completed_steps,ok,maximum_steps,absolute_tolerances)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: initial_energy_density,external_species(:),external_energy,duration,rtol,atol
    real(dp), intent(inout) :: species_density(:),temperature
    integer, intent(out) :: completed_steps
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_steps
    ! Species density tolerances followed by temperature tolerance; atol caps each entry.
    real(dp), optional, intent(in) :: absolute_tolerances(:)
    real(c_double), allocatable, target :: state(:),constraints(:),weights(:)
    real(dp) :: rates(size(species)),recovered,full_state(size(species)),final_density,final_temperature
    real(c_double) :: returned(1)
    integer(c_long) :: steps(1)
    integer(c_int), target :: slot
    integer(c_int) :: status,cleanup
    integer :: i,n,limit
    logical :: valid,accepted
    type(c_ptr) :: sunctx,memory
    type(N_Vector), pointer :: vector,limits,tolerances
    type(SUNMatrix), pointer :: matrix
    type(SUNLinearSolver), pointer :: solver
    ok=.false.;completed_steps=0;limit=100000
    if(present(maximum_steps)) limit=maximum_steps
    if(.not.all(ieee_is_finite([duration,rtol,atol]))) return
    if(duration<0.or.min(rtol,atol)<=0.or.limit<1) return
    if(present(absolute_tolerances)) then
      if(size(absolute_tolerances)/=size(species)+1) return
      if(.not.all(ieee_is_finite(absolute_tolerances))) return
      if(any(absolute_tolerances<=0._dp)) return
    endif
    call forced_reactor_rhs(species,reactions,initial_energy_density,0._dp,species_density,external_species, &
      external_energy,temperature,rates,recovered,valid)
    if(.not.valid) return
    if(duration==0) then
      ok=.true.;return
    endif
    final_density=sum(species_density)+duration*sum(external_species)
    if(.not.ieee_is_finite(final_density).or.final_density<=0) return
    if(size(species)==1) then
      full_state=final_density
      call forced_reactor_rhs(species,reactions,initial_energy_density,duration,full_state,external_species, &
        external_energy,recovered,rates,final_temperature,valid)
      if(.not.valid) return
      species_density=full_state;temperature=final_temperature;completed_steps=1;ok=.true.;return
    end if
    slot=0
    do i=1,size(registry)
      if(registry(i)%busy) cycle
      slot=int(i,c_int);exit
    enddo
    if(slot==0) return
    registry(slot)%busy=.true.;registry(slot)%species=species;registry(slot)%reactions=reactions
    registry(slot)%source=external_species;registry(slot)%energy=initial_energy_density
    registry(slot)%energy_source=external_energy;registry(slot)%guess=recovered
    registry(slot)%density=sum(species_density);registry(slot)%density_source=sum(external_species)
    registry(slot)%dependent=maxloc(species_density,dim=1)
    registry(slot)%independent=pack([(i,i=1,size(species))],[(i/=registry(slot)%dependent,i=1,size(species))])
    n=size(species)+1;allocate(state(n),constraints(n))
    state(:n-1)=species_density;state(n)=recovered;constraints=1
    sunctx=c_null_ptr;memory=c_null_ptr;nullify(vector,limits,tolerances,matrix,solver);accepted=.false.
    solve: block
      status=FSUNContext_Create(SUN_COMM_NULL,sunctx)
      if(status/=0.or..not.c_associated(sunctx)) exit solve
      vector=>FN_VMake_Serial(int(n,index_kind),state,sunctx)
      limits=>FN_VMake_Serial(int(n,index_kind),constraints,sunctx)
      if(.not.associated(vector).or..not.associated(limits)) exit solve
      matrix=>FSUNDenseMatrix(int(n,index_kind),int(n,index_kind),sunctx)
      if(.not.associated(matrix)) exit solve
      solver=>FSUNLinSol_Dense(vector,matrix,sunctx)
      if(.not.associated(solver)) exit solve
      memory=FCVodeCreate(CV_BDF,sunctx)
      if(.not.c_associated(memory)) exit solve
      status=FCVodeInit(memory,c_funloc(full_cvode_rhs),0._c_double,vector)
      if(status/=CV_SUCCESS) exit solve
      status=FCVodeSetUserData(memory,c_loc(slot))
      if(status/=CV_SUCCESS) exit solve
      if(present(absolute_tolerances)) then
        allocate(weights(n))
        weights=real(min(atol,absolute_tolerances),c_double)
        tolerances=>FN_VMake_Serial(int(n,index_kind),weights,sunctx)
        if(.not.associated(tolerances)) exit solve
        status=FCVodeSVtolerances(memory,real(rtol,c_double),tolerances)
      else
        status=FCVodeSStolerances(memory,real(rtol,c_double),real(atol,c_double))
      endif
      if(status/=CV_SUCCESS) exit solve
      ! A frozen SDC/MOL source can remove an initially absent species before
      ! chemistry produces it later in this interval. Such auxiliary forced
      ! trajectories need signed internal steps (as in native PelePhysics).
      ! Constraining every CVODE substep can stall at t=0 even when the final
      ! state is admissible. Keep constraints for nonnegative forcing; every
      ! endpoint still passes the same nonnegativity/EOS/rollback checks below.
      if(all(external_species>=0._dp)) then
        status=FCVodeSetConstraints(memory,limits)
        if(status/=CV_SUCCESS) exit solve
      endif
      status=FCVodeSetLinearSolver(memory,solver,matrix)
      if(status/=CV_SUCCESS) exit solve
      ! Pinned PelePhysics uses an independent-temperature approximate
      ! Jacobian although its RHS recovers temperature from affine energy.
      status=FCVodeSetJacFn(memory,c_funloc(full_cvode_jacobian))
      if(status/=CV_SUCCESS) exit solve
      ! Pinned CPU ReactorCvode::initCvode refreshes the approximate Jacobian
      ! after at most 100 steps; the CVODE default refreshes it earlier.
      status=FCVodeSetJacEvalFrequency(memory,100_c_long)
      if(status/=CV_SUCCESS) exit solve
      ! Match the pinned CPU reactor order and ordinary dense-output time policy.
      status=FCVodeSetMaxOrd(memory,2_c_int)
      if(status/=CV_SUCCESS) exit solve
      status=FCVodeSetMaxNumSteps(memory,int(limit,c_long))
      if(status/=CV_SUCCESS) exit solve
      status=FCVode(memory,real(duration,c_double),vector,returned,CV_NORMAL)
      if(status/=CV_SUCCESS.and.status/=CV_TSTOP_RETURN) exit solve
      if(abs(real(returned(1),dp)-duration)>100*epsilon(1._dp)*duration) exit solve
      status=FCVodeGetNumSteps(memory,steps)
      if(status/=CV_SUCCESS.or.steps(1)>int(limit,c_long)) exit solve
      ! Constraint-satisfying solves can return signed roundoff for absent
      ! trace species. Apply the existing machine-scale bound, then check mass
      ! independently of integration tolerances. No density projection is applied.
      if(any(state(:n-1) < -16*epsilon(1._dp)*final_density)) exit solve
      full_state=max(real(state(:n-1),dp),0._dp)
      if(abs(sum(full_state)-final_density)>5.e-12_dp*final_density) exit solve
      call forced_reactor_rhs(species,reactions,initial_energy_density,duration,full_state,external_species, &
        external_energy,real(state(n),dp),rates,recovered,valid)
      if(.not.valid) exit solve
      accepted=.true.
    end block solve
    if(c_associated(memory)) call FCVodeFree(memory)
    if(associated(solver)) then
      cleanup=FSUNLinSolFree(solver)
      if(cleanup/=0) accepted=.false.
    endif
    if(associated(matrix)) call FSUNMatDestroy(matrix)
    if(associated(tolerances)) call FN_VDestroy(tolerances)
    if(associated(limits)) call FN_VDestroy(limits)
    if(associated(vector)) call FN_VDestroy(vector)
    if(c_associated(sunctx)) then
      cleanup=FSUNContext_Free(sunctx)
      if(cleanup/=0) accepted=.false.
    endif
    deallocate(registry(slot)%species,registry(slot)%reactions,registry(slot)%source,registry(slot)%independent)
    registry(slot)%busy=.false.
    if(.not.accepted) return
    species_density=full_state;temperature=recovered;completed_steps=int(steps(1));ok=.true.
  end subroutine

  ! Full-coordinate temperature is an inverse-EOS seed, as in pinned
  ! PelePhysics fKernelSpec. The affine energy, not that coordinate, sets
  ! the rates and heat capacity at every internal CVODE evaluation.
  subroutine evaluate_forced_reactor_full_rhs(species,reactions,initial_energy,elapsed, &
      density_state,source,energy_source,temperature_seed,derivative,temperature_rate,temperature,ok)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: initial_energy,elapsed,density_state(:),source(:),energy_source,temperature_seed
    real(dp), intent(out) :: derivative(:),temperature_rate,temperature
    logical, intent(out) :: ok
    real(dp) :: rates(size(species)),recovered,cp,cv,h,e,entropy,capacity,thermal_source
    integer :: c
    logical :: valid
    ok=.false.;derivative=0;temperature_rate=0;temperature=0
    if(size(derivative)/=size(species)) return
    call forced_reactor_rhs(species,reactions,initial_energy,elapsed,density_state,source,energy_source, &
      temperature_seed,rates,recovered,valid,permit_signed_trials=.true.)
    if(.not.valid) return
    capacity=0;thermal_source=energy_source
    do c=1,size(species)
      call nasa7_mass_properties(species(c),recovered,cp,cv,h,e,entropy,valid)
      if(.not.valid) return
      capacity=capacity+density_state(c)*cv
      thermal_source=thermal_source-rates(c)*e
    enddo
    if(.not.ieee_is_finite(capacity).or.capacity<=0) return
    if(.not.ieee_is_finite(thermal_source)) return
    temperature_rate=thermal_source/capacity
    if(.not.ieee_is_finite(temperature_rate)) then
      temperature_rate=0;return
    endif
    derivative=rates;temperature=recovered;ok=.true.
  end subroutine evaluate_forced_reactor_full_rhs

  integer(c_int) function full_cvode_rhs(time,y,dy,user_data) result(status) bind(C)
    real(c_double), value :: time
    type(N_Vector) :: y,dy
    type(c_ptr), value :: user_data
    integer(c_int), pointer :: slot
    real(c_double), pointer :: state(:),rhs(:)
    real(dp) :: recovered,temperature_rate
    real(dp), allocatable :: rates(:)
    integer :: ns
    logical :: valid
    status=-1
    if(.not.c_associated(user_data)) return
    call c_f_pointer(user_data,slot)
    if(slot<1.or.slot>size(registry)) return
    if(.not.registry(slot)%busy) return
    state=>FN_VGetArrayPointer(y);rhs=>FN_VGetArrayPointer(dy)
    if(.not.associated(state).or..not.associated(rhs)) return
    ns=size(registry(slot)%species)
    if(size(state)/=ns+1.or.size(rhs)/=ns+1) return
    if(.not.all(ieee_is_finite(state))) return
    allocate(rates(ns));status=1
    call evaluate_forced_reactor_full_rhs(registry(slot)%species,registry(slot)%reactions, &
      registry(slot)%energy,real(time,dp),real(state(:ns),dp),registry(slot)%source, &
      registry(slot)%energy_source,real(state(ns+1),dp),rates,temperature_rate,recovered,valid)
    if(.not.valid) return
    rhs(:ns)=rates;rhs(ns+1)=temperature_rate;status=0
  end function full_cvode_rhs
  ! Only the approximate Jacobian uses this independent-temperature kernel.
  integer(c_int) function full_cvode_independent_temperature_rhs(time,y,dy,user_data) result(status) bind(C)
    real(c_double), value :: time
    type(N_Vector) :: y,dy
    type(c_ptr), value :: user_data
    integer(c_int), pointer :: slot
    real(c_double), pointer :: state(:),rhs(:)
    real(dp) :: recovered,cp,cv,h,e,entropy,capacity,thermal_source
    integer :: c,ns
    real(dp), allocatable :: full_state(:),rates(:)
    logical :: valid
    status=-1
    if(.not.c_associated(user_data)) return
    call c_f_pointer(user_data,slot)
    if(slot<1.or.slot>size(registry)) return
    if(.not.registry(slot)%busy) return
    state=>FN_VGetArrayPointer(y);rhs=>FN_VGetArrayPointer(dy)
    if(.not.associated(state).or..not.associated(rhs)) return
    allocate(full_state(size(registry(slot)%species)),rates(size(registry(slot)%species)))
    if(.not.all(ieee_is_finite(state))) return
    ! Newton trial coordinates may be negative before CVODE enforces its
    ! solution constraints. Evaluate their signed polynomial/thermodynamic state;
    ! the accepted endpoint still passes strict density/thermodynamic checks.
    full_state=real(state(:size(full_state)),dp)
    recovered=real(state(size(state)),dp)
    status=1
    if(.not.ieee_is_finite(recovered).or.recovered<=0) return
    block
      real(dp) :: concentrations(size(full_state)),production(size(full_state))
      do c=1,size(full_state)
        concentrations(c)=full_state(c)/registry(slot)%species(c)%molecular_weight
      enddo
      production=0
      if(size(registry(slot)%reactions)>0) then
        call elementary_production_rates_from_concentrations(registry(slot)%species,registry(slot)%reactions, &
          recovered,concentrations,production,valid,signed_concentrations=.true.)
        if(.not.valid) return
      endif
      do c=1,size(full_state)
        rates(c)=production(c)*registry(slot)%species(c)%molecular_weight+registry(slot)%source(c)
      enddo
    end block
    ns=size(full_state);rhs(:ns)=rates
    capacity=0;thermal_source=registry(slot)%energy_source
    do c=1,ns
      call nasa7_mass_properties(registry(slot)%species(c),recovered,cp,cv,h,e,entropy,valid)
      if(.not.valid) return
      capacity=capacity+full_state(c)*cv
      thermal_source=thermal_source-rates(c)*e
    enddo
    if(capacity<=0) return
    rhs(ns+1)=thermal_source/capacity
    status=0
  end function
  integer(c_int) function full_cvode_jacobian(time,y,dy,matrix,user_data,tmp1,tmp2,tmp3) &
      result(status) bind(C)
    real(c_double), value :: time
    type(N_Vector) :: y,dy,tmp1,tmp2,tmp3
    type(SUNMatrix) :: matrix
    type(c_ptr), value :: user_data
    integer(c_int), pointer :: slot
    real(c_double), pointer :: state(:),rhs(:),work(:),plus(:),minus(:),a(:,:)
    real(dp), allocatable :: conc(:),jac(:,:),energies(:),cvs(:)
    real(dp) :: cp,cv,h,e,entropy,capacity,delta
    integer :: ns,i,j
    logical :: valid
    status=-1
    if(.not.c_associated(user_data)) return
    call c_f_pointer(user_data,slot)
    if(slot<1.or.slot>size(registry)) return
    if(.not.registry(slot)%busy) return
    ns=size(registry(slot)%species)
    state=>FN_VGetArrayPointer(y);rhs=>FN_VGetArrayPointer(dy)
    work=>FN_VGetArrayPointer(tmp1);plus=>FN_VGetArrayPointer(tmp2);minus=>FN_VGetArrayPointer(tmp3)
    if(.not.associated(state).or..not.associated(rhs).or..not.associated(work)) return
    if(.not.associated(plus).or..not.associated(minus)) return
    if(size(state)/=ns+1.or..not.all(ieee_is_finite(state))) return
    a(1:ns+1,1:ns+1)=>FSUNDenseMatrix_Data(matrix)
    if(.not.associated(a)) return
    allocate(conc(ns),jac(ns,ns),energies(ns),cvs(ns))
    do i=1,ns
      conc(i)=state(i)/registry(slot)%species(i)%molecular_weight
      call nasa7_mass_properties(registry(slot)%species(i),state(ns+1),cp,cv,h,e,entropy,valid)
      if(.not.valid) return
      energies(i)=e;cvs(i)=cv
    enddo
    call elementary_production_jacobian(registry(slot)%species,registry(slot)%reactions, &
      state(ns+1),conc,jac,valid,signed_concentrations=.true.)
    if(.not.valid) return
    capacity=sum(state(:ns)*cvs)
    if(capacity<=0) return
    do j=1,ns
      do i=1,ns
        a(i,j)=jac(i,j)*registry(slot)%species(i)%molecular_weight/registry(slot)%species(j)%molecular_weight
      enddo
      a(ns+1,j)=(-sum(energies*a(:ns,j))-rhs(ns+1)*cvs(j))/capacity
    enddo
    ! Hybrid Jacobian: analytic density columns; centered temperature column.
    delta=max(1.e-4_dp,1.e-6_dp*state(ns+1))
    work=state;work(ns+1)=state(ns+1)+delta
    status=full_cvode_independent_temperature_rhs(time,tmp1,tmp2,user_data)
    if(status/=0) return
    work=state;work(ns+1)=state(ns+1)-delta
    status=full_cvode_independent_temperature_rhs(time,tmp1,tmp3,user_data)
    if(status/=0) return
    a(:,ns+1)=(plus-minus)/(2*delta)
    status=0
  end function full_cvode_jacobian
end module
