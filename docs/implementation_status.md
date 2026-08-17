# Implementation status

## Phase 0 — infrastructure

- [x] CMake and GNU Fortran builds
- [x] Debug bounds and floating-point checks
- [x] CTest and GitHub Actions
- [x] architecture, mapping, parity, and design-decision records

## Phase 1 — minimal Euler solver

- [x] uniform one-dimensional mesh
- [x] constant-`gamma` ideal-gas EOS
- [x] primitive/conserved conversion
- [x] outflow and periodic boundaries
- [x] Rusanov flux and conservative divergence
- [x] CFL control and SSPRK2
- [x] Sod exact-solution regression
- [x] entropy-wave convergence regression
- [x] Shu-Osher deterministic regression
- [x] symmetric planar Sedov-type strong-blast regression
- [ ] multidimensional isentropic-vortex regression

## Phase 2 — higher-order Godunov path

- [x] selectable PCM baseline
- [x] componentwise primitive PLM
- [x] minmod and MC limiters
- [x] selectable Rusanov/PeleC-style Riemann solver
- [x] characteristic projection and inverse mapping
- [x] one-dimensional PeleC characteristic tracing
- [x] dedicated time-centered conservative Godunov update
- [x] PeleC order-2 and order-4 limited slopes
- [x] pressure/velocity shock flattening
- [x] flattening unit gates for smooth, compressive, and expansive data
- [x] fourth-order-stencil smooth convergence gate
- [x] strong-shock positivity, symmetry, conservation, and signature gate
- [ ] general-EOS/internal-energy characteristic terms
- [ ] multispecies and passive-scalar tracing
- [ ] multidimensional transverse corrections
- [ ] PPM
- [ ] WENO

## Current verified strong-shock metrics

GNU Fortran Debug run, 800 cells, `t=0.02`:

| Metric | Result |
|---|---:|
| Completed steps | `1052` |
| Minimum density | `1.4352255e-1` |
| Maximum density | `5.0001945` |
| Minimum pressure | `1.0e-5` |
| Maximum pressure | `1.3654134e1` |
| Shock radius | `1.31875e-1` |
| Mass-balance error | `1.58e-14` |
| Momentum-balance error | `3.94e-33` |
| Energy-balance error | `3.29e-14` |
| Density/pressure symmetry error | `0` |
| Velocity antisymmetry error | `0` |

## Next implementation slice

Introduce a two-dimensional uniform Cartesian Euler scaffold and transverse Godunov corrections while retaining all one-dimensional baselines. Use a periodic isentropic vortex to verify multidimensional convergence before adding AMR or reacting-flow physics.
