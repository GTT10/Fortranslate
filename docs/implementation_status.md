# Implementation status

## Phase 0 — infrastructure

- [x] CMake and GNU Fortran builds
- [x] Debug bounds/FPE checks
- [x] CTest and GitHub Actions
- [x] architecture, mapping, parity, and design-decision records

## Phase 1 — minimal Euler solver

- [x] uniform one-dimensional mesh
- [x] constant-`gamma` ideal-gas EOS
- [x] primitive/conserved conversion
- [x] outflow and periodic boundaries
- [x] Rusanov flux
- [x] conservative flux divergence and CFL control
- [x] SSPRK2
- [x] Sod exact-solution regression
- [x] entropy-wave convergence regression
- [x] Shu-Osher deterministic regression
- [ ] Sedov regression
- [ ] isentropic-vortex regression

## Phase 2 — higher-order Godunov path

- [x] selectable PCM baseline
- [x] componentwise primitive PLM
- [x] minmod and MC limiters
- [x] positivity scaling and first-order face fallback
- [x] selectable Rusanov/PeleC-style Riemann solver
- [x] qualified single-species constant-`gamma` PeleC Riemann subset
- [x] characteristic-variable projection and inverse mapping
- [x] qualified one-dimensional PeleC characteristic tracing
- [x] dedicated time-centered conservative Godunov update
- [x] characteristic tracing unit tests
- [x] characteristic entropy-wave convergence gate
- [x] characteristic Sod exact-solution gate
- [x] characteristic Shu-Osher signature and conservation gate
- [ ] PeleC fourth-order slope option
- [ ] flattening
- [ ] general-EOS/internal-energy characteristic terms
- [ ] multispecies/passive-scalar tracing
- [ ] multidimensional transverse corrections
- [ ] PPM
- [ ] WENO

## Verified characteristic-PLM results

GNU Fortran 14.2 Debug build with bounds and floating-point traps:

| Test | Result |
|---|---:|
| Entropy wave, 40 cells, density L1 | `3.4864e-4` |
| Entropy wave, 80 cells, density L1 | `7.8837e-5` |
| Entropy wave, 160 cells, density L1 | `1.7751e-5` |
| Observed order, 40→80 | `2.1448` |
| Observed order, 80→160 | `2.1509` |
| Sod density L1 | `1.2178e-3` |
| Sod pressure L1 | `7.1822e-4` |
| Shu-Osher minimum density | `8.0000e-1` |
| Shu-Osher maximum density | `4.6137` |
| Shu-Osher interaction-window extrema | `19` |
| Shu-Osher mass-balance error | `8.17e-14` |
| Shu-Osher energy-balance error | `5.12e-13` |

GitHub Actions Debug and Release jobs pass the complete test suite.

## Next implementation slice

Add PeleC-compatible fourth-order limited slopes and a separately testable shock-flattening coefficient. Retain second-order MC slopes as the baseline, then add a Sedov regression before extending the characteristic state to general EOS or multispecies variables.
