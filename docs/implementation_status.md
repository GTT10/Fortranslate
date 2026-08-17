# Implementation status

## Phase 0 — project and test infrastructure

- [x] CMake build
- [x] GNU Fortran support
- [x] Debug bounds/FPE checks
- [x] CTest harness
- [x] GitHub Actions workflow
- [x] architecture, mapping, state, and parity documents

## Phase 1 — minimal compressible Euler solver

- [x] uniform one-dimensional mesh
- [x] constant-gamma ideal-gas EOS
- [x] primitive/conserved conversion
- [x] outflow ghost-cell boundaries
- [x] Rusanov flux
- [x] conservative finite-volume divergence
- [x] hydro CFL timestep
- [x] SSPRK2 time integration
- [x] Sod input and initialization
- [x] CSV output
- [x] exact Sod comparison
- [x] smooth-advection/entropy-wave regression
- [ ] Shu–Osher regression
- [ ] Sedov regression
- [ ] isentropic-vortex regression

## Phase 2 — higher-order Godunov path

- [x] retain selectable piecewise-constant baseline
- [x] componentwise primitive-variable PLM
- [x] minmod limiter
- [x] monotonized-central limiter
- [x] density/pressure face-state protection and first-order fallback
- [x] periodic ghost-cell boundaries
- [x] periodic slope wrapping
- [x] linear-reconstruction unit test
- [x] second-order smooth entropy-wave convergence test
- [x] tightened PLM Sod regression
- [ ] characteristic-variable projection
- [ ] PeleC-compatible characteristic tracing
- [ ] PeleC approximate Riemann solver
- [ ] multidimensional transverse corrections
- [ ] PPM
- [ ] WENO

## Verified PLM results

GNU Fortran 14.2 Debug build with bounds and floating-point traps:

| Test | Result |
|---|---:|
| Entropy wave, 40 cells, density L1 | `4.8297e-4` |
| Entropy wave, 80 cells, density L1 | `1.1659e-4` |
| Entropy wave, 160 cells, density L1 | `2.7383e-5` |
| Observed order, 40→80 | `2.0505` |
| Observed order, 80→160 | `2.0901` |
| Sod PLM density L1 | `1.8907e-3` |
| Sod PLM pressure L1 | `1.1981e-3` |

## Next implementation slice

Add a separately selectable PeleC-style approximate Riemann flux while retaining Rusanov as the robustness baseline. Validate equal-state consistency, shock-tube behavior, and a new Shu–Osher regression before characteristic tracing is attempted.
