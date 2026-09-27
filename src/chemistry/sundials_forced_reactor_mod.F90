module sundials_forced_reactor_mod
  use, intrinsic :: iso_c_binding
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use constant_volume_reactor_mod, only: forced_reactor_rhs
  use fsundials_core_mod
  use fcvode_mod
  use fnvector_serial_mod, only: FN_VMake_Serial
  use fsunmatrix_dense_mod, only: FSUNDenseMatrix
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
  public :: advance_forced_reactor_cvode
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
      status=FCVodeSetConstraints(memory,limits)
      if(status/=CV_SUCCESS) exit solve
      status=FCVodeSetLinearSolver(memory,solver,matrix)
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
    ! solution constraints. Extend the RHS as in the constant-volume backend;
    ! the accepted endpoint still passes strict density/thermodynamic checks.
    call reconstruct_density(registry(slot),real(time,dp),max(state,0._c_double),full_state)
    call forced_reactor_rhs(registry(slot)%species,registry(slot)%reactions,registry(slot)%energy, &
      real(time,dp),full_state,registry(slot)%source,registry(slot)%energy_source, &
      registry(slot)%guess,rates,recovered,valid)
    status=1
    if(.not.valid) return
    rhs=rates(registry(slot)%independent)
    registry(slot)%guess=recovered;status=0
  end function
  subroutine reconstruct_density(context,time,reduced,full_state)
    type(forcing_context), intent(in) :: context
    real(dp), intent(in) :: time
    real(c_double), intent(in) :: reduced(:)
    real(dp), intent(out) :: full_state(:)
    full_state(context%independent)=real(reduced,dp)
    full_state(context%dependent)=context%density+time*context%density_source-sum(real(reduced,dp))
  end subroutine
end module
