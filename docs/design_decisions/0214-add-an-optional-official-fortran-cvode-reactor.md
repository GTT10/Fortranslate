# 0214: Add an optional official-Fortran CVODE selected reactor

## Status

Accepted for `0.221.0`.

## Context

The `0.220.0` selected-mechanism application advances a configure-time
normalized bundle with the native adaptive backward-Euler reactor. That path
is deterministic and useful as a baseline, but it is not the Phase-7
production-stiff integration boundary. The pinned PelePhysics baseline carries
SUNDIALS commit `0eff39663606f2ff280c4059a947ed62ae38180a`, which identifies
itself as SUNDIALS `7.2.0` and provides an official Fortran 2003 module
interface when built with `BUILD_FORTRAN_MODULE_INTERFACE=ON`.

The reactor's thermodynamic state is not C interoperable: `nasa7_species` and
`elementary_reaction` contain ordinary Fortran data and allocatable members.
They must not cross a C ABI. The existing energy-constrained Jacobian is also
an `(N-1) x (N-1)` derivative, not the generated fixed-temperature `N x N`
kinetics Jacobian.

## Decision

Add `PELEF_ENABLE_SUNDIALS`, defaulting to `OFF`. When enabled, require exactly
SUNDIALS `7.2.0` and the official `fcvode_mod`, `fnvecserial_mod`,
`fsunmatrixdense_mod`, and `fsunlinsoldense_mod` CMake components and their
static imported targets. Shared-only installations are rejected so the
installed executable has no unshipped SUNDIALS runtime dependency. Compile
the adapter in a separate `pelef_selected_cvode_runtime` library so the
ordinary selected runtime, fixed applications, and mechanism probe remain
independent of SUNDIALS.

Keep one installed `pelef0d_selected` application and select its backend with
`integrator = "implicit"` or `integrator = "cvode"`. The default remains
`implicit`; a CVODE request in a build without the optional dependency fails
before CSV creation and never falls back silently. `maximum_steps` retains its
native implicit meaning, while `cvode_max_internal_steps` limits CVODE's
cumulative internal work across all output calls. Before each `FCVode` call,
the adapter supplies only the remaining budget to SUNDIALS.

Use CVODE BDF, the serial N_Vector, a dense SUNMatrix and dense linear solver,
scalar mass-fraction tolerances, nonnegative independent-species constraints,
and an explicit Jacobian callback. The CVODE state has `N-1` mass fractions.
At initialization, fix the largest initial mass-fraction species as the
dependent species and reconstruct it as
`Y_dep = 1 - sum(Y_independent)`. Every callback rejects nonfinite or invalid
trial states with a recoverable return; it does not publish them.

Evaluate the RHS through `reactor_rhs`. Generalize
`reactor_reduced_jacobian` to accept any dependent-species index and pass that
energy-constrained semi-analytic Jacobian directly to CVODE. Its kinetics
composition derivative is analytic, while its temperature derivative remains
a bounded centered finite difference. The generated mechanism rate kernel is
still used only for reported `wdot_*` columns.

Own SUNDIALS resources and private state in the adapter. Because the
noninteroperable mechanism context remains in module storage, allow exactly one
active reactor per process and declare the API serial, non-reentrant, and not
thread safe. Bind only the SUNDIALS callbacks to C. No species or reaction
derived type crosses the C boundary.

An advance commits caller `Y` and `T` only after CVODE returns the requested
time, the reconstructed composition is finite and admissible, fixed-energy
temperature recovery succeeds, and the total internal-step limit is not
exceeded. Any integration failure leaves caller state unchanged and marks the
private handle unusable until finalization. A scheduler fragment shorter than
`minimum_time_step` temporarily lowers CVODE's internal minimum to that exact
fragment, preserving the selected-reactor final-time contract.

## Consequences

The qualified optional path is a configure-time selected, serial, adiabatic,
constant-volume reactor for schema-1 mechanisms with 2--32 species and the
already supported elementary, third-body, and Troe forms. The two-species
fixture exercises a nonfinal dependent species, rollback, nested-context
rejection, repeat determinism, and a final fragment below the nominal minimum.
The pinned 10-species/29-reaction H2/O2 trajectory supplies the stiff,
conservation, and Cantera gates.

This decision does not add runtime YAML/JSON loading, multiple concurrent
reactors, OpenMP safety, sparse linear solvers, CVODES sensitivities,
selected-mechanism CFD/AMR/EB/MPI dispatch, detailed-fuel qualification, or a
performance claim. SUNDIALS 5/6, C-only, or shared-only installations are
intentionally rejected rather than adapted through a compatibility shim or
an incomplete runtime-install policy.
