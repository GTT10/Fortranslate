# PeleF architecture

## Current scope

The first implementation milestone is a serial one-dimensional finite-volume solver for the compressible Euler equations. It deliberately excludes AMR, chemistry, viscosity, particles, MPI, and accelerator support.

The executable currently performs this path:

```text
namelist input
    ↓
uniform mesh
    ↓
Sod initial condition
    ↓
primitive → conserved conversion
    ↓
CFL timestep
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
- `src/hydro`: state conversion, boundary conditions, Riemann fluxes, and finite-volume operators.
- `src/driver`: configuration, timestep selection, integration, and diagnostics.
- `src/problems`: problem-specific initial conditions.
- `src/io`: output without coupling physics kernels to a file format.
- `app`: executable assembly only.
- `tests`: unit and regression tests.
- `tools`: independent comparison and parity utilities.

Physics kernels operate on contiguous arrays and use no AMReX or PelePhysics data structures. The Phase-1 state array is organized as `state(variable, cell)`, with one ghost cell at each domain boundary.

## Numerical method in the first milestone

The solver advances

\[
\frac{\partial U}{\partial t} + \frac{\partial F(U)}{\partial x}=0
\]

with a cell-centered finite-volume discretization. Interface fluxes use the local Lax–Friedrichs/Rusanov formula

\[
F_{i+1/2}=\frac{F(U_L)+F(U_R)}{2}
-\frac{a_{\max}}{2}(U_R-U_L),
\]

where `a_max = max(|u_L| + c_L, |u_R| + c_R)`. Time integration uses two-stage SSPRK2. Reconstruction is piecewise constant in this milestone; PLM is a later, separately tested phase.

## Design constraints

1. The numerical code must remain independently testable without the application driver.
2. Invalid density or pressure is reported rather than silently repaired.
3. The finite-volume update must remain conservative.
4. Derived quantities such as pressure are computed through the EOS layer.
5. Future multidimensional, multispecies, AMR, and parallel layers must not require rewriting the basic EOS or Riemann APIs.
