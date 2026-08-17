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
- [x] outflow and periodic ghost-cell boundaries
- [x] Rusanov flux
- [x] conservative finite-volume divergence
- [x] hydro CFL timestep
- [x] SSPRK2 time integration
- [x] Sod input, initialization, CSV output, and exact comparison
- [x] smooth-advection/entropy-wave regression
- [x] Shu-Osher problem and deterministic regression
- [ ] Sedov regression
- [ ] isentropic-vortex regression

## Phase 2 — higher-order Godunov path

- [x] retain selectable piecewise-constant baseline
- [x] componentwise primitive-variable PLM
- [x] minmod and monotonized-central limiters
- [x] density/pressure face-state protection and first-order fallback
- [x] periodic slope wrapping
- [x] linear-reconstruction unit test
- [x] second-order smooth entropy-wave convergence tests
- [x] selectable Rusanov/PeleC-style Riemann dispatch
- [x] single-species constant-gamma reduction of PeleC approximate Riemann logic
- [x] equal-state, stationary-contact, Sod-interface, and invalid-selection unit gates
- [x] tightened PeleC-flux Sod regression
- [x] Shu-Osher shock-density-wave regression
- [ ] characteristic-variable projection
- [ ] PeleC-compatible characteristic tracing
- [ ] flattening
- [ ] multidimensional transverse corrections
- [ ] PPM
- [ ] WENO

## Verified results

GNU Fortran 14.2 Debug build with bounds and floating-point traps:

| Test | Result |
|---|---:|
| PeleC-flux entropy wave, 40 cells, density L1 | `4.6967e-4` |
| PeleC-flux entropy wave, 80 cells, density L1 | `1.1084e-4` |
| PeleC-flux entropy wave, 160 cells, density L1 | `2.5946e-5` |
| Observed order, 40→80 | `2.0831` |
| Observed order, 80→160 | `2.0949` |
| Sod PeleC-flux density L1 | `1.3678e-3` |
| Sod PeleC-flux pressure L1 | `8.0552e-4` |
| Shu-Osher minimum density | `8.0000e-1` |
| Shu-Osher maximum density | `4.6106` |
| Shu-Osher interaction-window extrema | `21` |
| Shu-Osher mass-balance error | `1.14e-13` |
| Shu-Osher energy-balance error | `1.25e-12` |

GitHub Actions Debug and Release jobs pass the complete unit and regression suite.

## Next implementation slice

Introduce characteristic-variable projection and a separately testable characteristic PLM path. Keep componentwise PLM, Rusanov, and the current PeleC-style solver available so changes in reconstruction can be isolated from changes in flux evaluation.
