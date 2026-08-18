# PeleF architecture

## Executable split

PeleF currently exposes two serial uniform-grid drivers that share the EOS, state conversion, limiter, characteristic tracing, and Riemann modules.

```text
pelef
  └─ one-dimensional PCM / PLM / PeleC-style Godunov paths

pelef2d
  └─ two-dimensional periodic CTU-style Euler path
```

The split keeps the verified one-dimensional implementation unchanged while multidimensional infrastructure is developed behind separate data layouts and regression gates.

## Shared state and physics

Both drivers advance

```text
(rho, rho*u, rho*v, rho*w, rho*E)
```

and use the same constant-`gamma` ideal-gas EOS. The primitive layout remains

```text
(rho, u, v, w, p).
```

The 1D state is stored as `state(variable, cell)` with ghost cells. The 2D state is stored as `state(variable, i, j)` over interior periodic cells; periodic neighbor indices are resolved explicitly by the CTU operator.

## One-dimensional path

```text
namelist
  ↓
Sod | Shu-Osher | Sedov initialization
  ↓
PCM | componentwise PLM | characteristic PeleC PLM
  ↓
Rusanov | qualified PeleC Riemann solver
  ↓
SSPRK2 or time-centered Godunov update
  ↓
physical-state, conservation, and problem-specific gates
```

The characteristic path supports order-2 or five-point order-4 slopes and optional shock flattening.

## Two-dimensional path

```text
2D namelist and uniform Cartesian mesh
  ↓
periodic isentropic-vortex initialization
  ↓
primitive conversion and limited x/y slopes
  ↓
normal characteristic tracing
  ├─ x: direct u-c, u, u+c tracing
  └─ y: rotate (u,v), trace with the x kernel, rotate back
  ↓
provisional x/y Riemann fluxes
  ↓
half-step transverse conservative corrections
  ├─ x-face states corrected by y-flux divergence
  └─ y-face states corrected by x-flux divergence
  ↓
final x/y Riemann fluxes
  ↓
one unsplit conservative update
  ↓
positivity and periodic conservation gates
```

`directional_flux_mod` owns the x/y momentum rotation, so neither Riemann solver contains duplicated y-direction algebra. `ctu_2d_mod` owns reconstruction, provisional fluxes, transverse correction, final fluxes, and the update.

## Transverse-correction safety

Each normal predictor state is physical before transverse correction. If the full half-step transverse update would make density or pressure invalid, the correction is scaled along the line from the base face state to the proposed corrected state. A bounded bisection determines the largest accepted factor `theta` in `[0,1]`.

The smooth isentropic-vortex regressions use `theta = 1` everywhere. The limiter exists as a defensive invariant for later stronger multidimensional cases.

## Dimensional reduction

A dedicated unit test repeats a one-dimensional periodic entropy wave across every y row, advances it with the 2D CTU operator, and compares it with the existing one-dimensional characteristic Godunov step at the same `dt`. The states agree to roundoff, verifying that transverse machinery vanishes for dimensionally reduced data.

## Design constraints

1. Existing 1D baselines remain unchanged and selectable.
2. Directional rotation is explicit and unit tested.
3. One numerical flux is shared by both cells adjacent to each face.
4. Periodic flux differences telescope, preserving all conserved integrals to roundoff.
5. Invalid density or pressure is rejected or corrected only by an explicitly measured positivity factor.
6. “PeleC-style” denotes a documented regular-grid subset, not full PeleC parity.
7. Future multispecies, transport, AMR, EB, and parallel layers must preserve the current EOS and flux-dispatch boundaries.

## Multispecies extension boundary

`multispecies_state_mod` owns the dynamic state layout and composition checks. `multispecies_flux_mod` and `reconstruction_multispecies_mod` add passive face transport without changing the existing five-variable Riemann APIs. The 2D hydro CTU routine may optionally return its provisional/final face states, fluxes, and positivity factors through `ctu_face_data_2d`; `ctu_multispecies_2d_mod` consumes that data so species and total-mass updates use identical face mass fluxes.
