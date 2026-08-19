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
- [x] molecular-weight database subset for H2, H, O, O2, OH, H2O, and N2
- [x] mixture molecular weight and gas constant
- [x] mixture `cp`, `cv`, `gamma`, enthalpy, and internal energy
- [x] ideal-gas pressure/density and frozen sound speed
- [x] bracketed Newton/bisection `e -> T` inversion
- [x] reference-value tests on both NASA7 coefficient intervals
- [ ] couple mixture thermodynamics into hydro state conversion and Riemann solvers
- [ ] composition-dependent hydro CFL and characteristic relations

## Phase 4 — zero-dimensional chemistry

- [x] synthetic first-order isomerization gate
- [x] arbitrary elementary reactant/product stoichiometry
- [x] reversible Arrhenius rates and NASA7 equilibrium constants
- [x] concentrations, progress rates, production rates, and `dY/dt`
- [x] normalized JSON-to-Fortran mechanism generator
- [x] committed generated-source cleanliness tests
- [x] adaptive explicit RK4 constant-volume reactor
- [x] stage-wise adiabatic `e -> T` coupling
- [x] seven-species, four-reaction H2/O2/N2 subset
- [x] third-body efficiencies
- [x] Lindemann reduced-pressure falloff
- [x] Troe broadening
- [x] duplicate reactions
- [x] analytic fixed-temperature production and mass-fraction Jacobians
- [x] reduced constant-energy reactor Jacobian
- [x] dense backward Euler/Newton integration with line search
- [x] adaptive step doubling and Richardson extrapolation
- [x] ten-species, 29-reaction H2/O2/Ar/N2 mechanism
- [x] structural mass, element, positivity, closure, and energy gates
- [x] live Cantera 3.2 trajectory and exact-state production-rate parity
- [ ] SRI and chemically activated reaction forms
- [ ] direct Cantera/CHEMKIN YAML parser
- [ ] sparse Jacobian and sparse linear solve
- [ ] CVODE/ARKODE integration
- [ ] detailed hydrocarbon mechanism
- [ ] chemistry coupling to the flow solver

## Verified `0.10.0` results

The ordinary suite contains 55 tests. Enabling the live Cantera reference adds the elementary and full-mechanism parity tests.

| Gate | Result |
|---|---:|
| Local Debug suite | `55/55` passed |
| Local Release suite | `55/55` passed |
| Generated elementary/full mechanism cleanliness | passed |
| Full H2/O2 output rows | `101` |
| Full H2/O2 accepted implicit steps | `5461` |
| Full H2/O2 final temperature | `2906.522001 K` |
| Full H2/O2 final pressure | `262444.241 Pa` |
| Maximum relative internal-energy error | `9.53e-12` |
| Maximum composition-closure error | `2.22e-16` |
| Maximum H inventory drift | `8.19e-16` |
| Maximum O inventory drift | `5.52e-16` |
| Maximum instantaneous mass-source residual | `4.68e-14` |
| Cantera-reference maximum temperature difference | `4.59e-3 K` |
| Cantera-reference maximum pressure difference | `3.72e-1 Pa` |
| Cantera-reference maximum species relative difference | `< 1.0e-5` |
| PeleF/Cantera final-temperature difference | `7.18e-5 K` |

The trajectory figures above use a pinned Cantera 3.2 constant-volume reference at the same 101 output times. The live CI gate also compares production rates at the exact PeleF `(T,rho,Y)` state so time-integration error is not confused with rate-kernel error.

## Next implementation slice

Couple NASA7 thermodynamics and the full chemistry source into a separate one-dimensional conservative-flow path. Begin with a first-order general-EOS Rusanov baseline and Strang splitting, then add general-EOS PLM/Godunov reconstruction and a non-uniform reacting-wave regression.
