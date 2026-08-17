# PeleF architecture

## Current scope

The current milestone is a serial one-dimensional finite-volume solver for the compressible Euler equations. It excludes AMR, chemistry, viscosity, particles, MPI, embedded boundaries, and accelerator support.

```text
namelist input
    ↓
problem selection and uniform mesh
    ↓
Sod or Shu-Osher initialization
    ↓
boundary fill
    ↓
reconstruction dispatch
    ├─ pcm: cell averages
    ├─ plm: componentwise primitive slopes
    └─ pelec_plm: primitive slopes + characteristic wave tracing
    ↓
Riemann dispatch
    ├─ Rusanov
    └─ PeleC-style ideal-gas subset
    ↓
conservative flux divergence
    ↓
time update
    ├─ SSPRK2 for pcm/plm
    └─ one time-centered Godunov update for pelec_plm
    ↓
physical-state gate, CSV output, independent regressions
```

## Layering

- `src/core`: precision, constants, state indices, and mesh utilities.
- `src/physics`: EOS functions independent of discretization.
- `src/hydro`: state conversion, boundaries, limiting, reconstruction, characteristic tracing, Riemann solvers, and finite-volume operators.
- `src/driver`: configuration, timestep selection, integration dispatch, and diagnostics.
- `src/problems`: problem-specific initial conditions.
- `src/io`: output isolated from physics kernels.
- `app`: executable assembly and problem selection only.
- `tests`: algebraic, convergence, exact-solution, and deterministic regression gates.
- `tools`: independent Python comparison utilities.

The state array remains `state(variable, cell)` with one ghost cell on each side. Physics kernels use contiguous Fortran arrays and no AMReX or PelePhysics data structures.

## Replaceable numerical components

- `reconstruction = "pcm" | "plm" | "pelec_plm"`;
- `limiter = "minmod" | "mc"`;
- `boundary_condition = "outflow" | "periodic"`;
- `riemann_solver = "rusanov" | "pelec"`.

Unknown selections fail explicitly. A selected PeleC-style component is never silently replaced with Rusanov or first order.

`pelec_plm` differs structurally from the method-of-lines paths because its interface states already include `dt/dx` characteristic tracing. The driver therefore uses one conservative Godunov update for that path instead of applying SSPRK2 to an already time-centered flux.

## Design constraints

1. Every numerical component must remain independently testable.
2. Invalid density or pressure is rejected rather than silently repaired.
3. Flux-divergence updates remain conservative.
4. First-order, componentwise-PLM, and Rusanov baselines remain available for differential diagnosis.
5. “PeleC-style” always names a documented and tested subset, not whole-solver parity.
6. General EOS, multispecies, multidimensional, AMR, and parallel extensions must not require replacing the existing dispatch boundaries.
