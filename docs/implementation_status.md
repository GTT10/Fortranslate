# Implementation status

## Phase 0 — infrastructure

- [x] CMake and GNU Fortran builds
- [x] Debug bounds and floating-point checks
- [x] CTest and GitHub Actions
- [x] architecture, mapping, parity, and design-decision records

## Phase 1 — minimal Euler solver

- [x] one-dimensional uniform mesh
- [x] constant-`gamma` ideal-gas EOS
- [x] primitive/conserved conversion
- [x] outflow and periodic 1D boundaries
- [x] Rusanov flux and conservative divergence
- [x] CFL control and SSPRK2
- [x] Sod exact-solution regression
- [x] entropy-wave convergence regression
- [x] Shu-Osher deterministic regression
- [x] symmetric planar Sedov-type strong-blast regression

## Phase 2 — higher-order Godunov path

- [x] selectable PCM baseline
- [x] componentwise primitive PLM
- [x] minmod and MC limiters
- [x] selectable Rusanov/PeleC-style Riemann solver
- [x] characteristic projection and inverse mapping
- [x] one-dimensional PeleC characteristic tracing
- [x] PeleC order-2 and order-4 limited slopes
- [x] pressure/velocity shock flattening
- [x] dedicated time-centered conservative 1D Godunov update
- [x] two-dimensional uniform Cartesian state and mesh scaffold
- [x] x/y directional flux rotation
- [x] normal characteristic tracing in both directions
- [x] provisional multidimensional fluxes
- [x] conservative half-step transverse flux corrections
- [x] positivity scaling for transverse corrections
- [x] unsplit conservative 2D update
- [x] dimensional-reduction parity with the 1D solver
- [x] periodic isentropic-vortex analytical regression
- [x] multidimensional second-order convergence gate
- [ ] physical 2D wall, inflow, and outflow boundaries
- [ ] general-EOS/internal-energy characteristic terms
- [ ] PPM
- [ ] WENO
- [ ] 3D double-transverse corrections

## Verified 2D results

GNU Fortran Debug build with bounds and floating-point traps:

| Test | Result |
|---|---:|
| Vortex density L1, `24 x 24` | `2.5862222e-3` |
| Vortex density L1, `48 x 48` | `5.3334804e-4` |
| Vortex density L1, `96 x 96` | `1.1010295e-4` |
| Observed order, 24→48 | `2.2777` |
| Observed order, 48→96 | `2.2762` |
| 48² density L1 without transverse correction | `1.1493797e-3` |
| Maximum periodic conservation error | `5.68e-14` |
| Minimum transverse positivity factor | `1.0` |

The complete Debug and Release suites contain 31 tests and pass without failures.

## Phase 3 — multispecies advection

- [ ] configurable species count and state indexing
- [ ] conserved species densities `rho*Y_k`
- [ ] species fluxes through all Riemann paths
- [ ] positivity and normalization handling
- [ ] one-dimensional MultiSpecSod regression
- [ ] two-dimensional passive/species dimensional-reduction gate

## Next implementation slice

Introduce multispecies conserved-state infrastructure and non-reacting species advection. Preserve the current five-variable Euler kernel as the zero-species specialization, then add a MultiSpecSod parity case before NASA thermodynamics or chemistry is attempted.
