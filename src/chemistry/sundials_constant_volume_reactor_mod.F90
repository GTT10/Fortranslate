module sundials_constant_volume_reactor_mod
  use, intrinsic :: iso_c_binding, only: &
    c_associated, c_double, c_f_pointer, c_funloc, c_int, c_int32_t, &
    c_int64_t, c_loc, c_long, c_null_ptr, c_ptr
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use mixture_thermo_mod, only: &
    temperature_from_internal_energy, valid_mixture_composition
  use elementary_kinetics_mod, only: elementary_reaction
  use constant_volume_reactor_mod, only: &
    reactor_rhs, reactor_reduced_jacobian
  use fsundials_core_mod, only: &
    N_Vector, SUNLinearSolver, SUNMatrix, SUN_COMM_NULL, &
    FN_VDestroy, FN_VGetArrayPointer, FSUNContext_Create, &
    FSUNContext_Free, FSUNLinSolFree, FSUNMatDestroy
  use fcvode_mod, only: &
    CV_BDF, CV_NORMAL, CV_SUCCESS, FCVode, FCVodeCreate, FCVodeFree, &
    FCVodeGetNumErrTestFails, FCVodeGetNumJacEvals, &
    FCVodeGetNumNonlinSolvConvFails, FCVodeGetNumNonlinSolvIters, &
    FCVodeGetNumRhsEvals, FCVodeGetNumSteps, FCVodeInit, &
    FCVodeSetConstraints, FCVodeSetInitStep, FCVodeSetJacFn, &
    FCVodeSetLinearSolver, FCVodeSetMaxNumSteps, FCVodeSetMaxStep, &
    FCVodeSetMinStep, FCVodeSetUserData, FCVodeSStolerances
  use fnvector_serial_mod, only: FN_VMake_Serial
  use fsunmatrix_dense_mod, only: FSUNDenseMatrix, FSUNDenseMatrix_Data
  use fsunlinsol_dense_mod, only: FSUNLinSol_Dense
  implicit none
  private

#ifdef SUNDIALS_INT32_T
  integer, parameter :: sunindex_kind = c_int32_t
#else
  integer, parameter :: sunindex_kind = c_int64_t
#endif
  integer, parameter, public :: cvode_reactor_context_capacity = 64

  type, public :: cvode_reactor_context
    private
    integer :: slot = 0
    integer(int64) :: generation = 0_int64
  end type cvode_reactor_context

  type, public :: cvode_reactor_statistics
    integer(int64) :: internal_steps = 0_int64
    integer(int64) :: rhs_evaluations = 0_int64
    integer(int64) :: jacobian_evaluations = 0_int64
    integer(int64) :: nonlinear_iterations = 0_int64
    integer(int64) :: nonlinear_convergence_failures = 0_int64
    integer(int64) :: error_test_failures = 0_int64
  end type cvode_reactor_statistics

  type, bind(C) :: cvode_callback_token
    integer(c_int) :: slot = 0_c_int
    integer(c_int64_t) :: generation = 0_c_int64_t
  end type cvode_callback_token

  type :: cvode_context_storage
    logical :: active = .false.
    logical :: failed = .false.
    integer :: dependent_species = 0
    integer :: maximum_steps = 0
    integer(int64) :: generation = 0_int64
    integer, allocatable :: independent_species(:)
    real(dp) :: density = 0.0_dp
    real(dp) :: target_energy = 0.0_dp
    real(dp) :: temperature_guess = 0.0_dp
    real(dp) :: time = 0.0_dp
    real(dp) :: minimum_step = 0.0_dp
    type(nasa7_species), allocatable :: species(:)
    type(elementary_reaction), allocatable :: reactions(:)
    real(dp), allocatable :: full_state(:)
    real(dp), allocatable :: full_rhs(:)
    real(dp), allocatable :: jacobian(:, :)
    real(c_double), pointer :: reduced_state(:) => null()
    real(c_double), pointer :: constraints(:) => null()
    type(N_Vector), pointer :: solution_vector => null()
    type(N_Vector), pointer :: constraint_vector => null()
    type(SUNMatrix), pointer :: system_matrix => null()
    type(SUNLinearSolver), pointer :: linear_solver => null()
    type(c_ptr) :: cvode_memory = c_null_ptr
    type(c_ptr) :: sundials_context = c_null_ptr
  end type cvode_context_storage

  type(cvode_context_storage), target, save :: &
    context_registry(cvode_reactor_context_capacity)
  type(cvode_callback_token), target, save :: &
    callback_tokens(cvode_reactor_context_capacity)
  integer(int64), save :: next_context_generation = 1_int64

  public :: initialize_constant_volume_cvode
  public :: advance_constant_volume_cvode
  public :: get_constant_volume_cvode_statistics
  public :: finalize_constant_volume_cvode

contains

  subroutine initialize_constant_volume_cvode( &
      context, species, reactions, density, target_internal_energy, initial_time, &
      initial_time_step, minimum_time_step, maximum_time_step, &
      relative_tolerance, absolute_tolerance, maximum_steps, &
      mass_fractions, temperature, ok, message)
    type(cvode_reactor_context), intent(inout) :: context
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: density, target_internal_energy, initial_time
    real(dp), intent(in) :: initial_time_step, minimum_time_step
    real(dp), intent(in) :: maximum_time_step
    real(dp), intent(in) :: relative_tolerance, absolute_tolerance
    integer, intent(in) :: maximum_steps
    real(dp), intent(in) :: mass_fractions(:), temperature
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    real(dp), allocatable :: initial_rhs(:)
    real(dp) :: recovered_temperature
    type(cvode_context_storage), pointer :: storage
    integer(c_int) :: status
    integer(sunindex_kind) :: reduced_size
    integer :: i, independent_index, slot

    ok = .false.
    message = ""
    if (context%slot /= 0 .or. context%generation /= 0_int64) then
      message = "CVODE reactor handle is already assigned"
      return
    end if
    if (size(species) < 2 .or. size(mass_fractions) /= size(species)) then
      message = "CVODE reactor state size is invalid"
      return
    end if
    if (.not. valid_mixture_composition(species, mass_fractions)) then
      message = "CVODE reactor initial composition is invalid"
      return
    end if
    if (density <= 0.0_dp .or. relative_tolerance <= 0.0_dp .or. &
        absolute_tolerance <= 0.0_dp .or. initial_time_step <= 0.0_dp .or. &
        minimum_time_step <= 0.0_dp .or. maximum_time_step < minimum_time_step .or. &
        initial_time_step < minimum_time_step .or. &
        initial_time_step > maximum_time_step .or. maximum_steps <= 0 .or. &
        .not. all(ieee_is_finite([ &
          density, target_internal_energy, initial_time, initial_time_step, &
          minimum_time_step, maximum_time_step, relative_tolerance, &
          absolute_tolerance, temperature]))) then
      message = "CVODE reactor scalar configuration is invalid"
      return
    end if

    allocate(initial_rhs(size(species)))
    call reactor_rhs( &
      species, reactions, density, target_internal_energy, mass_fractions, &
      temperature, initial_rhs, recovered_temperature, ok)
    if (.not. ok) then
      message = "CVODE reactor could not recover the initial thermodynamic state"
      return
    end if
    ok = .false.

    slot = find_available_context_slot()
    if (slot == 0) then
      message = "CVODE reactor context capacity is exhausted"
      return
    end if
    storage => context_registry(slot)
    storage%active = .true.
    storage%failed = .false.
    storage%generation = next_context_generation
    if (next_context_generation == huge(next_context_generation)) then
      next_context_generation = 1_int64
    else
      next_context_generation = next_context_generation + 1_int64
    end if
    callback_tokens(slot)%slot = int(slot, c_int)
    callback_tokens(slot)%generation = int(storage%generation, c_int64_t)
    storage%species = species
    storage%reactions = reactions
    storage%density = density
    storage%target_energy = target_internal_energy
    storage%temperature_guess = recovered_temperature
    storage%time = initial_time
    storage%minimum_step = minimum_time_step
    storage%maximum_steps = maximum_steps
    storage%dependent_species = maxloc(mass_fractions, dim=1)
    allocate(storage%independent_species(size(species) - 1))
    independent_index = 0
    do i = 1, size(species)
      if (i == storage%dependent_species) cycle
      independent_index = independent_index + 1
      storage%independent_species(independent_index) = i
    end do
    allocate(storage%full_state(size(species)), storage%full_rhs(size(species)))
    allocate(storage%jacobian(size(species) - 1, size(species) - 1))
    allocate(storage%reduced_state(size(species) - 1))
    allocate(storage%constraints(size(species) - 1))
    do i = 1, size(storage%independent_species)
      storage%reduced_state(i) = real( &
        mass_fractions(storage%independent_species(i)), c_double)
    end do
    storage%constraints = 1.0_c_double

    status = FSUNContext_Create(SUN_COMM_NULL, storage%sundials_context)
    if (status /= 0_c_int .or. &
        .not. c_associated(storage%sundials_context)) then
      call fail_initialization(slot, "FSUNContext_Create", status, message)
      return
    end if
    reduced_size = int(size(storage%reduced_state), sunindex_kind)
    storage%solution_vector => FN_VMake_Serial( &
      reduced_size, storage%reduced_state, storage%sundials_context)
    if (.not. associated(storage%solution_vector)) then
      call fail_initialization( &
        slot, "FN_VMake_Serial(solution)", -1_c_int, message)
      return
    end if
    storage%constraint_vector => FN_VMake_Serial( &
      reduced_size, storage%constraints, storage%sundials_context)
    if (.not. associated(storage%constraint_vector)) then
      call fail_initialization( &
        slot, "FN_VMake_Serial(constraints)", -1_c_int, message)
      return
    end if
    storage%system_matrix => FSUNDenseMatrix( &
      reduced_size, reduced_size, storage%sundials_context)
    if (.not. associated(storage%system_matrix)) then
      call fail_initialization(slot, "FSUNDenseMatrix", -1_c_int, message)
      return
    end if
    storage%linear_solver => FSUNLinSol_Dense( &
      storage%solution_vector, storage%system_matrix, storage%sundials_context)
    if (.not. associated(storage%linear_solver)) then
      call fail_initialization(slot, "FSUNLinSol_Dense", -1_c_int, message)
      return
    end if
    storage%cvode_memory = FCVodeCreate(CV_BDF, storage%sundials_context)
    if (.not. c_associated(storage%cvode_memory)) then
      call fail_initialization(slot, "FCVodeCreate", -1_c_int, message)
      return
    end if

    status = FCVodeInit( &
      storage%cvode_memory, c_funloc(pelef_cvode_rhs_callback), &
      real(initial_time, c_double), storage%solution_vector)
    if (status /= CV_SUCCESS) then
      call fail_initialization(slot, "FCVodeInit", status, message)
      return
    end if
    status = FCVodeSetUserData( &
      storage%cvode_memory, c_loc(callback_tokens(slot)))
    if (status /= CV_SUCCESS) then
      call fail_initialization(slot, "FCVodeSetUserData", status, message)
      return
    end if
    status = FCVodeSStolerances( &
      storage%cvode_memory, real(relative_tolerance, c_double), &
      real(absolute_tolerance, c_double))
    if (status /= CV_SUCCESS) then
      call fail_initialization(slot, "FCVodeSStolerances", status, message)
      return
    end if
    status = FCVodeSetConstraints( &
      storage%cvode_memory, storage%constraint_vector)
    if (status /= CV_SUCCESS) then
      call fail_initialization(slot, "FCVodeSetConstraints", status, message)
      return
    end if
    status = FCVodeSetLinearSolver( &
      storage%cvode_memory, storage%linear_solver, storage%system_matrix)
    if (status /= CV_SUCCESS) then
      call fail_initialization(slot, "FCVodeSetLinearSolver", status, message)
      return
    end if
    status = FCVodeSetJacFn( &
      storage%cvode_memory, c_funloc(pelef_cvode_jacobian_callback))
    if (status /= CV_SUCCESS) then
      call fail_initialization(slot, "FCVodeSetJacFn", status, message)
      return
    end if
    status = FCVodeSetInitStep( &
      storage%cvode_memory, real(initial_time_step, c_double))
    if (status /= CV_SUCCESS) then
      call fail_initialization(slot, "FCVodeSetInitStep", status, message)
      return
    end if
    status = FCVodeSetMinStep( &
      storage%cvode_memory, real(minimum_time_step, c_double))
    if (status /= CV_SUCCESS) then
      call fail_initialization(slot, "FCVodeSetMinStep", status, message)
      return
    end if
    status = FCVodeSetMaxStep( &
      storage%cvode_memory, real(maximum_time_step, c_double))
    if (status /= CV_SUCCESS) then
      call fail_initialization(slot, "FCVodeSetMaxStep", status, message)
      return
    end if
    status = FCVodeSetMaxNumSteps( &
      storage%cvode_memory, int(maximum_steps, c_long))
    if (status /= CV_SUCCESS) then
      call fail_initialization(slot, "FCVodeSetMaxNumSteps", status, message)
      return
    end if

    context%slot = slot
    context%generation = storage%generation
    message = ""
    ok = .true.
  end subroutine initialize_constant_volume_cvode

  subroutine advance_constant_volume_cvode( &
      context, target_time, mass_fractions, temperature, ok, message)
    type(cvode_reactor_context), intent(in) :: context
    real(dp), intent(in) :: target_time
    real(dp), intent(inout) :: mass_fractions(:), temperature
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    type(cvode_reactor_statistics) :: statistics
    type(cvode_context_storage), pointer :: storage
    real(dp), allocatable :: candidate_state(:)
    real(dp) :: candidate_temperature, time_tolerance
    real(c_double) :: returned_time(1)
    integer(c_int) :: status
    integer(c_long) :: completed_steps(1), remaining_steps

    ok = .false.
    message = ""
    call resolve_context_handle(context, storage, ok, message)
    if (.not. ok) return
    ok = .false.
    if (storage%failed) then
      message = "The CVODE reactor context is failed and must be finalized"
      return
    end if
    if (size(mass_fractions) /= size(storage%species)) then
      message = "CVODE reactor output state has the wrong size"
      return
    end if
    if (.not. ieee_is_finite(target_time) .or. target_time <= storage%time) then
      message = "CVODE reactor target time must advance monotonically"
      return
    end if

    status = FCVodeGetNumSteps(storage%cvode_memory, completed_steps)
    if (status /= CV_SUCCESS) then
      call fail_advance(storage, "FCVodeGetNumSteps", status, message)
      return
    end if
    remaining_steps = int(storage%maximum_steps, c_long) - completed_steps(1)
    if (remaining_steps <= 0_c_long) then
      storage%failed = .true.
      message = "CVODE reactor cumulative step limit reached"
      return
    end if
    status = FCVodeSetMaxNumSteps(storage%cvode_memory, remaining_steps)
    if (status /= CV_SUCCESS) then
      call fail_advance(storage, "FCVodeSetMaxNumSteps", status, message)
      return
    end if

    status = FCVodeSetMinStep( &
      storage%cvode_memory, real(min( &
        storage%minimum_step, target_time - storage%time), c_double))
    if (status /= CV_SUCCESS) then
      call fail_advance(storage, "FCVodeSetMinStep", status, message)
      return
    end if
    returned_time = real(storage%time, c_double)
    status = FCVode( &
      storage%cvode_memory, real(target_time, c_double), &
      storage%solution_vector, &
      returned_time, CV_NORMAL)
    if (status /= CV_SUCCESS) then
      call fail_advance(storage, "FCVode", status, message)
      return
    end if
    time_tolerance = 100.0_dp * epsilon(1.0_dp) * &
      max(1.0_dp, abs(target_time))
    if (abs(real(returned_time(1), dp) - target_time) > time_tolerance) then
      storage%failed = .true.
      message = "CVODE did not return the requested output time"
      return
    end if

    allocate(candidate_state(size(storage%species)))
    call reconstruct_full_state( &
      storage, storage%reduced_state, candidate_state, ok)
    if (.not. ok) then
      storage%failed = .true.
      message = "CVODE returned an invalid mass-fraction state"
      return
    end if
    call temperature_from_internal_energy( &
      storage%species, candidate_state, storage%target_energy, &
      storage%temperature_guess, candidate_temperature, ok)
    if (.not. ok) then
      storage%failed = .true.
      message = "CVODE returned an invalid thermodynamic state"
      return
    end if
    call get_constant_volume_cvode_statistics( &
      context, statistics, ok, message)
    if (.not. ok) then
      storage%failed = .true.
      return
    end if
    if (statistics%internal_steps > int(storage%maximum_steps, int64)) then
      storage%failed = .true.
      message = "CVODE reactor step limit reached"
      ok = .false.
      return
    end if

    mass_fractions = candidate_state
    temperature = candidate_temperature
    storage%temperature_guess = candidate_temperature
    storage%time = target_time
    message = ""
    ok = .true.
  end subroutine advance_constant_volume_cvode

  subroutine get_constant_volume_cvode_statistics( &
      context, statistics, ok, message)
    type(cvode_reactor_context), intent(in) :: context
    type(cvode_reactor_statistics), intent(out) :: statistics
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    integer(c_long) :: value(1)
    integer(c_int) :: status
    type(cvode_context_storage), pointer :: storage

    statistics = cvode_reactor_statistics()
    ok = .false.
    message = ""
    call resolve_context_handle(context, storage, ok, message)
    if (.not. ok) return
    ok = .false.

    status = FCVodeGetNumSteps(storage%cvode_memory, value)
    if (status /= CV_SUCCESS) then
      call format_status("FCVodeGetNumSteps", status, message)
      return
    end if
    statistics%internal_steps = int(value(1), int64)
    status = FCVodeGetNumRhsEvals(storage%cvode_memory, value)
    if (status /= CV_SUCCESS) then
      call format_status("FCVodeGetNumRhsEvals", status, message)
      return
    end if
    statistics%rhs_evaluations = int(value(1), int64)
    status = FCVodeGetNumJacEvals(storage%cvode_memory, value)
    if (status /= CV_SUCCESS) then
      call format_status("FCVodeGetNumJacEvals", status, message)
      return
    end if
    statistics%jacobian_evaluations = int(value(1), int64)
    status = FCVodeGetNumNonlinSolvIters(storage%cvode_memory, value)
    if (status /= CV_SUCCESS) then
      call format_status("FCVodeGetNumNonlinSolvIters", status, message)
      return
    end if
    statistics%nonlinear_iterations = int(value(1), int64)
    status = FCVodeGetNumNonlinSolvConvFails(storage%cvode_memory, value)
    if (status /= CV_SUCCESS) then
      call format_status( &
        "FCVodeGetNumNonlinSolvConvFails", status, message)
      return
    end if
    statistics%nonlinear_convergence_failures = int(value(1), int64)
    status = FCVodeGetNumErrTestFails(storage%cvode_memory, value)
    if (status /= CV_SUCCESS) then
      call format_status("FCVodeGetNumErrTestFails", status, message)
      return
    end if
    statistics%error_test_failures = int(value(1), int64)
    ok = .true.
  end subroutine get_constant_volume_cvode_statistics

  subroutine finalize_constant_volume_cvode(context, ok, message)
    type(cvode_reactor_context), intent(inout) :: context
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    integer(c_int) :: status
    type(cvode_context_storage), pointer :: storage
    integer :: slot
    character(len=64) :: cleanup_stage

    if (context%slot == 0 .and. context%generation == 0_int64) then
      ok = .true.
      message = ""
      return
    end if
    call resolve_context_handle(context, storage, ok, message)
    if (.not. ok) then
      context%slot = 0
      context%generation = 0_int64
      return
    end if
    slot = context%slot
    call release_context(slot, status, cleanup_stage)
    context%slot = 0
    context%generation = 0_int64
    ok = status == 0_c_int
    if (ok) then
      message = ""
    else
      call format_status( &
        trim(cleanup_stage) // " during CVODE reactor cleanup", &
        status, message)
    end if
  end subroutine finalize_constant_volume_cvode

  integer(c_int) function pelef_cvode_rhs_callback( &
      time, state_vector, rhs_vector, user_data) result(status) &
      bind(C, name="pelef_cvode_rhs_callback")
    real(c_double), value :: time
    type(N_Vector) :: state_vector, rhs_vector
    type(c_ptr), value :: user_data

    real(c_double), pointer :: reduced_state(:), reduced_rhs(:)
    real(dp) :: recovered_temperature
    type(cvode_context_storage), pointer :: context
    logical :: state_ok
    integer :: i

    status = -1_c_int
    call resolve_callback_context(user_data, context, state_ok)
    if (.not. state_ok) return
    if (context%failed) return
    reduced_state => FN_VGetArrayPointer(state_vector)
    reduced_rhs => FN_VGetArrayPointer(rhs_vector)
    if (.not. associated(reduced_state) .or. &
        .not. associated(reduced_rhs)) return
    if (size(reduced_state) /= size(context%independent_species) .or. &
        size(reduced_rhs) /= size(context%independent_species)) return
    call reconstruct_full_state( &
      context, reduced_state, context%full_state, state_ok)
    if (.not. state_ok) then
      status = 1_c_int
      return
    end if
    call reactor_rhs( &
      context%species, context%reactions, context%density, &
      context%target_energy, context%full_state, context%temperature_guess, &
      context%full_rhs, recovered_temperature, state_ok)
    if (.not. state_ok) then
      status = 1_c_int
      return
    end if
    do i = 1, size(context%independent_species)
      reduced_rhs(i) = real( &
        context%full_rhs(context%independent_species(i)), c_double)
    end do
    context%temperature_guess = recovered_temperature
    if (.not. ieee_is_finite(real(time, dp))) then
      status = -1_c_int
    else
      status = 0_c_int
    end if
  end function pelef_cvode_rhs_callback

  integer(c_int) function pelef_cvode_jacobian_callback( &
      time, state_vector, rhs_vector, jacobian_matrix, user_data, &
      temporary_1, temporary_2, temporary_3) result(status) &
      bind(C, name="pelef_cvode_jacobian_callback")
    real(c_double), value :: time
    type(N_Vector) :: state_vector, rhs_vector
    type(SUNMatrix) :: jacobian_matrix
    type(c_ptr), value :: user_data
    type(N_Vector) :: temporary_1, temporary_2, temporary_3

    real(c_double), pointer :: reduced_state(:), matrix_values(:, :)
    real(dp) :: recovered_temperature
    type(cvode_context_storage), pointer :: context
    logical :: state_ok
    integer :: reduced_size

    status = -1_c_int
    call resolve_callback_context(user_data, context, state_ok)
    if (.not. state_ok) return
    if (context%failed) return
    reduced_state => FN_VGetArrayPointer(state_vector)
    if (.not. associated(reduced_state)) return
    reduced_size = size(context%independent_species)
    if (size(reduced_state) /= reduced_size) return
    call reconstruct_full_state( &
      context, reduced_state, context%full_state, state_ok)
    if (.not. state_ok) then
      status = 1_c_int
      return
    end if
    call reactor_reduced_jacobian( &
      context%species, context%reactions, context%density, &
      context%target_energy, context%full_state, context%temperature_guess, &
      context%jacobian, recovered_temperature, state_ok, &
      context%dependent_species)
    if (.not. state_ok) then
      status = 1_c_int
      return
    end if
    matrix_values(1:reduced_size, 1:reduced_size) => &
      FSUNDenseMatrix_Data(jacobian_matrix)
    if (.not. associated(matrix_values)) return
    matrix_values = real(context%jacobian, c_double)
    context%temperature_guess = recovered_temperature
    if (.not. ieee_is_finite(real(time, dp))) then
      status = -1_c_int
    else
      status = 0_c_int
    end if
    if (associated(FN_VGetArrayPointer(rhs_vector))) continue
    if (associated(FN_VGetArrayPointer(temporary_1))) continue
    if (associated(FN_VGetArrayPointer(temporary_2))) continue
    if (associated(FN_VGetArrayPointer(temporary_3))) continue
  end function pelef_cvode_jacobian_callback

  subroutine reconstruct_full_state(context, reduced_state, full_state, ok)
    type(cvode_context_storage), intent(in) :: context
    real(c_double), intent(in) :: reduced_state(:)
    real(dp), intent(out) :: full_state(:)
    logical, intent(out) :: ok

    real(dp), parameter :: tolerance = 5.0e-12_dp
    integer :: i

    full_state = 0.0_dp
    ok = size(full_state) == size(reduced_state) + 1 .and. &
      size(reduced_state) == size(context%independent_species)
    if (.not. ok) return
    if (any(.not. ieee_is_finite(real(reduced_state, dp)))) then
      ok = .false.
      return
    end if
    do i = 1, size(reduced_state)
      full_state(context%independent_species(i)) = real(reduced_state(i), dp)
    end do
    full_state(context%dependent_species) = &
      1.0_dp - sum(real(reduced_state, dp))
    if (any(full_state < -tolerance) .or. &
        any(full_state > 1.0_dp + tolerance)) then
      ok = .false.
      return
    end if
    ok = valid_mixture_composition(context%species, full_state)
  end subroutine reconstruct_full_state

  integer function find_available_context_slot() result(slot)
    integer :: candidate

    slot = 0
    do candidate = 1, cvode_reactor_context_capacity
      if (.not. context_registry(candidate)%active .and. &
          .not. has_allocated_context(context_registry(candidate))) then
        slot = candidate
        return
      end if
    end do
  end function find_available_context_slot

  subroutine resolve_context_handle(context, storage, ok, message)
    type(cvode_reactor_context), intent(in) :: context
    type(cvode_context_storage), pointer, intent(out) :: storage
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    nullify(storage)
    ok = .false.
    message = ""
    if (context%slot == 0 .and. context%generation == 0_int64) then
      message = "CVODE reactor handle is inactive"
      return
    end if
    if (context%slot < 1 .or. &
        context%slot > cvode_reactor_context_capacity .or. &
        context%generation == 0_int64) then
      message = "CVODE reactor handle is invalid"
      return
    end if
    storage => context_registry(context%slot)
    if (.not. storage%active .or. &
        storage%generation /= context%generation) then
      nullify(storage)
      message = "CVODE reactor handle is stale"
      return
    end if
    if (.not. c_associated(storage%cvode_memory)) then
      nullify(storage)
      message = "CVODE reactor context is incomplete"
      return
    end if
    ok = .true.
  end subroutine resolve_context_handle

  subroutine resolve_callback_context(user_data, context, ok)
    type(c_ptr), intent(in) :: user_data
    type(cvode_context_storage), pointer, intent(out) :: context
    logical, intent(out) :: ok

    type(cvode_callback_token), pointer :: token
    integer :: slot

    nullify(context)
    nullify(token)
    ok = .false.
    if (.not. c_associated(user_data)) return
    call c_f_pointer(user_data, token)
    if (.not. associated(token)) return
    slot = int(token%slot)
    if (slot < 1 .or. slot > cvode_reactor_context_capacity) return
    if (.not. context_registry(slot)%active) return
    if (int(token%generation, int64) /= &
        context_registry(slot)%generation) return
    context => context_registry(slot)
    ok = .true.
  end subroutine resolve_callback_context

  logical function has_allocated_context(context) result(allocated_context)
    type(cvode_context_storage), intent(in) :: context

    allocated_context = allocated(context%species) .or. &
      allocated(context%reactions) .or. &
      allocated(context%independent_species) .or. &
      allocated(context%full_state) .or. allocated(context%full_rhs) .or. &
      allocated(context%jacobian) .or. &
      associated(context%reduced_state) .or. &
      associated(context%constraints) .or. &
      associated(context%solution_vector) .or. &
      associated(context%constraint_vector) .or. &
      associated(context%system_matrix) .or. &
      associated(context%linear_solver) .or. &
      c_associated(context%cvode_memory) .or. &
      c_associated(context%sundials_context)
  end function has_allocated_context

  subroutine fail_initialization(slot, stage, status, message)
    integer, intent(in) :: slot
    character(len=*), intent(in) :: stage
    integer(c_int), intent(in) :: status
    character(len=*), intent(out) :: message

    integer(c_int) :: cleanup_status
    character(len=64) :: cleanup_stage
    character(len=len(message)) :: cleanup_message

    call format_status(stage, status, message)
    call release_context(slot, cleanup_status, cleanup_stage)
    if (cleanup_status /= 0_c_int) then
      call format_status( &
        trim(cleanup_stage) // " during failed initialization cleanup", &
        cleanup_status, cleanup_message)
      message = trim(message) // "; " // trim(cleanup_message)
    end if
  end subroutine fail_initialization

  subroutine fail_advance(context, stage, status, message)
    type(cvode_context_storage), intent(inout) :: context
    character(len=*), intent(in) :: stage
    integer(c_int), intent(in) :: status
    character(len=*), intent(out) :: message

    context%failed = .true.
    call format_status(stage, status, message)
  end subroutine fail_advance

  subroutine format_status(stage, status, message)
    character(len=*), intent(in) :: stage
    integer(c_int), intent(in) :: status
    character(len=*), intent(out) :: message

    write(message, '(a," failed with status ",i0)') trim(stage), status
  end subroutine format_status

  subroutine release_context(slot, status, failed_stage)
    integer, intent(in) :: slot
    integer(c_int), intent(out) :: status
    character(len=*), intent(out), optional :: failed_stage

    type(cvode_context_storage), pointer :: context
    integer(c_int) :: free_status

    status = 0_c_int
    if (present(failed_stage)) failed_stage = ""
    if (slot < 1 .or. slot > cvode_reactor_context_capacity) then
      status = -1_c_int
      if (present(failed_stage)) failed_stage = "release_context"
      return
    end if
    context => context_registry(slot)
    context%active = .false.
    context%failed = .false.
    if (c_associated(context%cvode_memory)) then
      call FCVodeFree(context%cvode_memory)
    end if
    context%cvode_memory = c_null_ptr
    callback_tokens(slot)%slot = 0_c_int
    callback_tokens(slot)%generation = 0_c_int64_t
    if (associated(context%linear_solver)) then
      free_status = FSUNLinSolFree(context%linear_solver)
      if (status == 0_c_int .and. free_status /= 0_c_int) then
        status = free_status
        if (present(failed_stage)) failed_stage = "FSUNLinSolFree"
      end if
      nullify(context%linear_solver)
    end if
    if (associated(context%system_matrix)) then
      call FSUNMatDestroy(context%system_matrix)
      nullify(context%system_matrix)
    end if
    if (associated(context%constraint_vector)) then
      call FN_VDestroy(context%constraint_vector)
      nullify(context%constraint_vector)
    end if
    if (associated(context%solution_vector)) then
      call FN_VDestroy(context%solution_vector)
      nullify(context%solution_vector)
    end if
    if (c_associated(context%sundials_context)) then
      free_status = FSUNContext_Free(context%sundials_context)
      if (status == 0_c_int .and. free_status /= 0_c_int) then
        status = free_status
        if (present(failed_stage)) failed_stage = "FSUNContext_Free"
      end if
    end if
    context%sundials_context = c_null_ptr
    if (allocated(context%species)) deallocate(context%species)
    if (allocated(context%reactions)) deallocate(context%reactions)
    if (allocated(context%independent_species)) then
      deallocate(context%independent_species)
    end if
    if (allocated(context%full_state)) deallocate(context%full_state)
    if (allocated(context%full_rhs)) deallocate(context%full_rhs)
    if (allocated(context%jacobian)) deallocate(context%jacobian)
    if (associated(context%reduced_state)) then
      deallocate(context%reduced_state)
      nullify(context%reduced_state)
    end if
    if (associated(context%constraints)) then
      deallocate(context%constraints)
      nullify(context%constraints)
    end if
    context%dependent_species = 0
    context%maximum_steps = 0
    context%density = 0.0_dp
    context%target_energy = 0.0_dp
    context%temperature_guess = 0.0_dp
    context%time = 0.0_dp
    context%minimum_step = 0.0_dp
  end subroutine release_context

end module sundials_constant_volume_reactor_mod
