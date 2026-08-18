# Implementation status

## Phase 0 — infrastructure

- [x] CMake and GNU Fortran builds
- [x] Debug bounds and floating-point traps
- [x] CTest and GitHub Actions
- [x] architecture, mapping, parity, and design-decision records

## Phase 1 — uniform-grid Euler core

- [x] one-dimensional constant-`gamma` Euler solver
- [x] Rusanov and qualified PeleC-style Riemann solvers
- [x] componentwise and characteristic PLM
- [x] order-2/order-4 limited slopes and shock flattening
- [x] Sod, Shu-Osher, and planar Sedov-type regressions
- [x] two-dimensional periodic CTU-style scaffold
- [x] directional flux rotation and transverse corrections
- [x] isentropic-vortex convergence and dimensional reduction

## Phase 2 — passive multispecies transport

- [x] runtime species layout with conserved `rho*Y_k`
- [x] positivity and species-closure validation
- [x] species flux closure against the shared mass flux
- [x] 1D characteristic and 2D CTU species transport
- [x] MultiSpecSod and periodic species-wave regressions

## Phase 3 — composition-dependent thermodynamics

- [x] NASA7 species record and polynomial evaluation
- [x] molecular-weight database subset for H2, O2, H2O, and N2
- [x] mixture molecular weight and gas constant
- [x] mixture `cp`, `cv`, `gamma`, enthalpy, and internal energy
- [x] ideal-gas pressure/density and frozen sound speed
- [x] bracketed Newton/bisection `e -> T` inversion
- [x] reference-value tests on both NASA7 coefficient intervals
- [ ] couple mixture thermodynamics into hydro state conversion and Riemann solvers
- [ ] composition-dependent hydro CFL and characteristic relations

## Phase 4 — zero-dimensional reaction scaffold

- [x] Arrhenius first-order mass-conserving reaction object
- [x] isothermal RK4 integration and analytical solution gate
- [x] adiabatic constant-volume stage-wise `e -> T` coupling
- [x] energy, composition, monotonicity, and application-level CSV gates
- [ ] general reaction stoichiometry and production-rate kernel
- [ ] reversible rates and equilibrium constants
- [ ] mechanism parser/code generator
- [ ] analytic or generated Jacobian
- [ ] stiff integrator and Cantera parity
- [ ] detailed combustion mechanism

## Verified `0.8.0` results

GNU Fortran 14.2 Debug build with bounds and floating-point traps:

| Gate | Result |
|---|---:|
| Complete Debug suite | `44/44` passed |
| NASA7 H2 `cp` at 300 K | `14311.7571456761 J/(kg K)` |
| NASA7 O2 `cp` at 1500 K | `1143.02005736187 J/(kg K)` |
| Air-subset `Wmix` | `28.850334 kg/kmol` |
| Air-subset `gamma` at 1200 K | `1.32225210530124` |
| Air-subset frozen sound speed at 1200 K | `676.222203367036 m/s` |
| 0D final temperature | `1839.99999999720 K` |
| 0D final reactant fraction | `2.13e-110` |
| Maximum 0D relative energy error | `9.38e-12` |
| Maximum 0D composition closure error | `1.11e-16` |

## Next implementation slice

Add a mechanism representation and source-code generator for arbitrary elementary Arrhenius reactions, then validate a real constant-volume H2/O2 reactor against Cantera before coupling chemistry to the multidimensional flow solver.
