# PeleF architecture

## Current scope

The current milestone is a serial one-dimensional finite-volume solver for the compressible Euler equations. It deliberately excludes AMR, chemistry, viscosity, particles, MPI, and accelerator support.

The executable performs this path:

```text
namelist input
    ↓
problem selection and uniform mesh
    ↓
Sod or Shu-Osher initialization
    ↓
boundary fill
    ↓
PCM or PLM reconstruction
    ↓
Riemann-solver dispatch
    ├─ Rusanov
    └─ PeleC-style ideal-gas subset
    ↓
finite-volume flux divergence
    ↓
SSPRK2 update and physical-state gate
    ↓
CSV output and independent regression tools
```

## Layering

- `src/core`: precision, constants, state indices, and mesh utilities.
- `src/physics`: equation-of-state implementations independent of discretization.
- `src/hydro`: state conversion, boundary conditions, slope limiting, reconstruction, Riemann solvers, dispatch, and finite-volume operators.
- `src/driver`: configuration, timestep selection, integration, and diagnostics.
- `src/problems`: problem-specific initial conditions.
- `src/io`: output without coupling physics kernels to a file format.
- `app`: executable assembly and problem selection only.
- `tests`: unit and regression tests.
- `tools`: independent comparison and parity utilities.

Physics kernels operate on contiguous arrays and use no AMReX or PelePhysics data structures. The current state array is organized as `state(variable, cell)`, with one ghost cell at each domain boundary.

## Replaceable numerical components

The finite-volume operator receives names for the reconstruction, limiter, boundary condition, and Riemann solver. Dispatch is explicit:

- `reconstruction = "pcm" | "plm"`;
- `limiter = "minmod" | "mc"`;
- `boundary_condition = "outflow" | "periodic"`;
- `riemann_solver = "rusanov" | "pelec"`.

Unknown values fail configuration or flux evaluation. The code does not silently replace a failed selected solver with Rusanov; Rusanov is an explicit robustness baseline.

## Boundary behavior

The boundary layer owns ghost filling independently of reconstruction:

- `outflow`: copy the nearest interior cell into each ghost cell and suppress boundary-adjacent PLM slopes;
- `periodic`: wrap ghost states and the corresponding PLM slopes.

Wrapping slopes as well as states is required to retain second-order convergence across the periodic interface.

## Design constraints

1. Numerical kernels remain independently testable without the application driver.
2. Invalid density or pressure is reported rather than silently repaired.
3. The finite-volume update remains conservative.
4. Derived quantities such as pressure are computed through the EOS layer.
5. First-order and Rusanov paths remain available after higher-order components are added.
6. A subsystem is described as PeleC-compatible only for the tested subset explicitly documented.
7. Future characteristic reconstruction, multidimensional, multispecies, AMR, and parallel layers must not require rewriting the current EOS or flux-dispatch APIs.
