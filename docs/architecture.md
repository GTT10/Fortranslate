# PeleF architecture

## Current scope

PeleF currently advances the one-dimensional compressible Euler equations on a uniform Cartesian mesh. The state layout remains `state(variable, cell)` with one conserved-state ghost cell on each side. Reconstruction builds a wider temporary primitive array when a method needs a larger stencil.

```text
namelist input
    ↓
problem selection: Sod | Shu-Osher | Sedov
    ↓
uniform mesh and conserved state
    ↓
boundary fill
    ↓
reconstruction dispatch
    ├─ pcm
    ├─ plm: componentwise primitive reconstruction + SSPRK2
    └─ pelec_plm
         ├─ order-2 or PeleC order-4 limited slopes
         ├─ optional pressure/velocity shock flattening
         ├─ characteristic projection
         └─ u-c, u, u+c time tracing
    ↓
Riemann dispatch: Rusanov | qualified PeleC subset
    ↓
shared conservative face fluxes
    ↓
time update and physical-state gate
```

## Replaceable numerical components

The application passes explicit selections through the driver to the finite-volume operator:

- `reconstruction = "pcm" | "plm" | "pelec_plm"`;
- `limiter = "minmod" | "mc"`;
- `plm_order = 2 | 4`;
- `use_flattening = .false. | .true.`;
- `boundary_condition = "outflow" | "periodic"`;
- `riemann_solver = "rusanov" | "pelec"`.

The new controls only affect `pelec_plm`. Existing inputs omit them and therefore retain the verified order-2, no-flattening behavior.

## Stencil ownership

`reconstruction_pelec_plm_mod` owns the extended primitive stencil:

- order-2 slopes require cells `i-1:i+1`;
- order-4 slopes require `i-2:i+2`;
- flattening requires pressure and normal velocity through `i-3:i+3`.

For periodic boundaries, these temporary primitive and slope values wrap. For outflow boundaries, exterior primitive values are constant extensions and boundary-adjacent slopes are suppressed. The conserved-state storage and existing boundary module therefore remain unchanged.

## Time integration

`pcm` and componentwise `plm` form a method-of-lines spatial operator and use SSPRK2. `pelec_plm` constructs time-centered characteristic face states with `dt/dx`, so it uses a single conservative Godunov update. Mixing SSPRK2 with already traced states is intentionally rejected.

## Design constraints

1. Each algorithmic component has a direct unit or regression gate.
2. Density and pressure failures are reported; no silent solver substitution is allowed.
3. Rusanov and lower-order reconstruction remain selectable diagnostic baselines.
4. One face flux is shared by neighboring cells, preserving finite-volume conservation.
5. “PeleC-style” always denotes an explicitly documented subset, not whole-code parity.
6. Future general-EOS, multispecies, multidimensional, AMR, and parallel layers must preserve the current component boundaries.
