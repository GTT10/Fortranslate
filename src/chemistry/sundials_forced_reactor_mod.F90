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
    real(dp) :: rates(size(species)),recovered
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
    slot=0
    do i=1,size(registry)
      if(registry(i)%busy) cycle
      slot=int(i,c_int);exit
    enddo
    if(slot==0) return
    registry(slot)%busy=.true.;registry(slot)%species=species;registry(slot)%reactions=reactions
    registry(slot)%source=external_species;registry(slot)%energy=initial_energy_density
    registry(slot)%energy_source=external_energy;registry(slot)%guess=recovered
    n=size(species);allocate(state(n),constraints(n));state=species_density;constraints=1
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
      ! CVODE forms a numerical Jacobian of the complete forced density RHS.
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
      call forced_reactor_rhs(species,reactions,initial_energy_density,duration,real(state,dp),external_species, &
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
    deallocate(registry(slot)%species,registry(slot)%reactions,registry(slot)%source)
    registry(slot)%busy=.false.
    if(.not.accepted) return
    species_density=real(state,dp);temperature=recovered;completed_steps=int(steps(1));ok=.true.
  end subroutine

  integer(c_int) function forced_cvode_rhs(time,y,dy,user_data) result(status) bind(C)
    real(c_double), value :: time
    type(N_Vector) :: y,dy
    type(c_ptr), value :: user_data
    integer(c_int), pointer :: slot
    real(c_double), pointer :: state(:),rhs(:)
    real(dp) :: recovered
    logical :: valid
    status=-1
    if(.not.c_associated(user_data)) return
    call c_f_pointer(user_data,slot)
    if(slot<1.or.slot>size(registry)) return
    if(.not.registry(slot)%busy) return
    state=>FN_VGetArrayPointer(y);rhs=>FN_VGetArrayPointer(dy)
    if(.not.associated(state).or..not.associated(rhs)) return
    call forced_reactor_rhs(registry(slot)%species,registry(slot)%reactions,registry(slot)%energy, &
      real(time,dp),real(state,dp),registry(slot)%source,registry(slot)%energy_source, &
      registry(slot)%guess,rhs,recovered,valid)
    status=1
    if(.not.valid) return
    registry(slot)%guess=recovered;status=0
  end function
end module
