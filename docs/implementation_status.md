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
- [x] committed generated-source cleanliness test
- [x] adaptive explicit RK4 constant-volume reactor
- [x] stage-wise adiabatic `e -> T` coupling
- [x] seven-species, four-reaction H2/O2/N2 subset
- [x] structural mass, element, positivity, closure, and energy gates
- [x] live Cantera 3.2 trajectory and production-rate parity
- [ ] third-body efficiencies
- [ ] pressure falloff and Troe/SRI forms
- [ ] direct Cantera/CHEMKIN YAML parser
- [ ] analytic or generated Jacobian
- [ ] stiff integration through CVODE or an equivalent solver
- [ ] complete H2/O2 mechanism
- [ ] detailed hydrocarbon mechanism
- [ ] chemistry coupling to the flow solver

## Verified `0.9.0` results

GNU Fortran builds and the live Cantera 3.2 reference gate:

| Gate | Result |
|---|---:|
| Complete Debug suite | `50/50` passed |
| Complete Release suite | `50/50` passed |
| Generated mechanism cleanliness | passed |
| H2/O2 accepted adaptive steps | `751` |
| H2/O2 output rows | `101` |
| Initial temperature | `1200 K` |
| Final temperature | `1332.56148597839 K` |
| Final pressure | `112518.160472301 Pa` |
| Maximum relative internal-energy error | `9.92e-12` |
| Maximum composition-closure error | `2.22e-16` |
| Maximum H inventory drift | `4.51e-17` |
| Maximum O inventory drift | `1.91e-17` |
| Maximum instantaneous mass-source residual | `7.11e-13` |
| Cantera maximum temperature difference | `1.61e-6 K` |
| Cantera maximum pressure difference | `1.36e-4 Pa` |
| Cantera maximum species mass-fraction difference | `1.70e-11` |
| Cantera maximum production-rate difference | `3.55e-12 kmol/(m^3 s)` |
| PeleF/Cantera final-temperature difference | `3.69e-9 K` |

The production-rate maximum occurs in an almost cancelled OH net source of approximately `2.5e-8 kmol/(m^3 s)`; the absolute comparison floor is therefore `5e-12 kmol/(m^3 s)` while the relative tolerance remains active for non-cancelled rates.

## Next implementation slice

Add third-body and falloff/Troe reaction forms, generate an analytic or semi-generated Jacobian, and connect a stiff integrator. Then expand from the four-reaction subset to a complete small H2/O2 mechanism before attempting chemistry-flow coupling.
