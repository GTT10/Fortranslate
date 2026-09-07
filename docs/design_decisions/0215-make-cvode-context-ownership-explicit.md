# 0215: Make CVODE context ownership explicit

## Status

Accepted for `0.222.0`.

## Context

Milestone `0.221.0` introduced the optional official-Fortran SUNDIALS 7.2.0
CVODE backend, but its mechanism data, work arrays, solver resources, clock,
statistics, and failure state were module-global. It therefore rejected a
second live reactor even when a caller needed only sequentially interleaved
cell or case integration. Moving the existing Fortran derived types directly
through CVODE user data is not valid: `nasa7_species` and
`elementary_reaction` contain allocatable, non-C-interoperable components.

The ownership boundary also needs to survive ordinary Fortran intrinsic
assignment. A shallow copy of raw `c_ptr` or pointer members would permit
double destruction. Reusing a released registry slot without a generation
check would let an old handle destroy an unrelated new reactor.

## Decision

Make every public operation context-first:

```fortran
type(cvode_reactor_context) :: context

call initialize_constant_volume_cvode(context, ...)
call advance_constant_volume_cvode(context, target_time, y, t, ok, message)
call get_constant_volume_cvode_statistics(context, stats, ok, message)
call finalize_constant_volume_cvode(context, ok, message)
```

Keep `cvode_reactor_context` opaque. It contains only a private integer slot
and 64-bit generation. A fixed module-private registry provides 64 live slots;
each slot owns the complete Fortran mechanism state, reduced and full work
arrays, SUNDIALS context, N_Vectors, dense matrix, dense linear solver, CVODE
memory, clock, temperature guess, cumulative internal-step budget, and failure
flag. Publish the capacity as `cvode_reactor_context_capacity` so callers can
bound storage and tests can exercise exhaustion deterministically.

Associate each slot with a stable module-lifetime `BIND(C)` token containing
only C-interoperable slot and generation integers. After `FCVodeInit`, install
the token address with `FCVodeSetUserData`. RHS and Jacobian callbacks convert
only that token, validate the range, active flag, and generation, and then
resolve the registry entry. No mechanism, reaction, allocatable array, or
Fortran context object is passed through the C ABI.

Initialization reserves an empty slot only after input and initial
thermodynamic validation. Its result remains false until every SUNDIALS setup
call succeeds; capacity or later setup failure releases partial resources and
leaves the requested handle empty. Each advance resolves exactly one handle.
Its clock, remaining internal-step allowance, statistics, rollback, and failed
quarantine are independent of all peers.

Finalization first marks the slot inactive, then frees CVODE memory, linear
solver, matrix, vectors, SUNContext, and Fortran arrays in dependency order.
It clears the callback token only after CVODE is gone and finally clears the
public handle. Repeating finalization on that cleared handle succeeds.
Intrinsic copies alias one live slot rather than copying resources; the first
valid finalization releases it, after which any surviving copy is rejected as
stale and cleared. A stale copy cannot release a new context that reuses the
same slot because its generation differs.

If a vendor destroy routine returns a nonzero status, continue the ordered
cleanup, invalidate the handle, and identify the first failing routine in the
diagnostic. Retrying a possibly partial vendor destruction is not a supported
operation. Require the SUNDIALS precision header to declare double precision
only because the official Fortran callback ABI used here is `c_double`.

## Consequences

Up to 64 reactors can coexist and their operations can be interleaved in any
sequential order. Failure or finalization of one does not contaminate a peer.
The numerical method, tolerances, dependent-species rule, cumulative budget,
CSV contract, static SUNDIALS dependency, and native-default behavior remain
unchanged.

The acceptance gate compares two different-density, different-temperature
states with different dependent species. Their standalone and interleaved
final states, temperatures, and all six solver counters must be exactly equal.
It also forces one context to fail while the other advances, finalizes in
arbitrary and repeated order, reuses a slot behind a stale copied handle, and
fills and releases the entire registry.

This is serial multi-context ownership, not simultaneous execution. Registry
allocation, resolution, callbacks, and teardown are not thread safe or
reentrant. Sparse solvers, performance/scaling, runtime mechanism loading,
selected-mechanism CFD/AMR/EB/MPI dispatch, mechanisms above 32 species,
detailed fuels, and physical validation remain outside this decision.
