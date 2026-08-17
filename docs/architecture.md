# PeleF architecture

## Current scope

The current milestone is a serial one-dimensional finite-volume solver for the compressible Euler equations. It deliberately excludes AMR, chemistry, viscosity, particles, MPI, and accelerator support.

The executable performs this path:

```text
namelist input
    ↓
uniform mesh
    ↓
Sod initial condition
    ↓
boundary fill
    ↓
PCM or PLM reconstruction
    ↓
Rusanov interface fluxes
    ↓
finite-volume divergence
    ↓
SSPRK2 update
    ↓
CSV output and exact-solution comparison
```

## Layering

- `src/core`: precision, constants, state indices, and mesh utilities.
- `src/physics`: equation-of-state implementations independent of discretization.
- `src/hydro`: state conversion, boundary conditions, slope limiting, reconstruction, Riemann fluxes, and finite-volume operators.
- `src/driver`: configuration, timestep selection, integration, and diagnostics.
- `src/problems`: problem-specific initial conditions.
- `src/io`: output without coupling physics kernels to a file format.
- `app`: executable assembly only.
- `tests`: unit and regression tests.
- `tools`: independent comparison and parity utilities.

Physics kernels operate on contiguous arrays and use no AMReX or PelePhysics data structures. The current state array is organized as `state(variable, cell)`, with one ghost cell at each domain boundary.

## Spatial discretization

The solver advances

\[
\frac{\partial U}{\partial t} + \frac{\partial F(U)}{\partial x}=0
\]

with a cell-centered finite-volume discretization. Two reconstruction paths are retained:

1. `pcm`: piecewise-constant cell states, kept as the first-order regression baseline.
2. `plm`: componentwise primitive-variable slopes with minmod or monotonized-central limiting.

PLM face states are converted back to conserved variables before entering the Riemann solver. Density and pressure slopes are scaled if a face extrapolation would approach their configured floors; failed face conversion falls back to the corresponding cell-centered conserved state.

Interface fluxes use the local Lax–Friedrichs/Rusanov formula

\[
F_{i+1/2}=\frac{F(U_L)+F(U_R)}{2}
-\frac{a_{\max}}{2}(U_R-U_L),
\]

where `a_max = max(|u_L| + c_L, |u_R| + c_R)`. Time integration uses two-stage SSPRK2.

## Boundary behavior

The boundary layer owns ghost filling independently of reconstruction:

- `outflow`: copy the nearest interior cell into each ghost cell and suppress boundary-adjacent PLM slopes;
- `periodic`: wrap ghost states and wrap the corresponding PLM slopes.

Wrapping slopes as well as states is required to retain second-order convergence across the periodic interface.

## Design constraints

1. Numerical kernels must remain independently testable without the application driver.
2. Invalid density or pressure is reported rather than silently repaired.
3. The finite-volume update must remain conservative.
4. Derived quantities such as pressure are computed through the EOS layer.
5. The first-order path remains available after higher-order methods are added.
6. Future characteristic reconstruction, multidimensional, multispecies, AMR, and parallel layers must not require rewriting the EOS or Riemann APIs.
7. Componentwise PLM is treated as a verified scaffold, not as full PeleC PLM parity.
