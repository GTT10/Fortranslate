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
- [x] qualified 1D conserved/primitive conversion using NASA7 mixture thermodynamics
- [x] composition-dependent 1D CFL and Rusanov signal speeds
- [x] frozen-composition characteristic relations
- [ ] general-EOS PeleC-style Riemann parity
- [ ] multidimensional general-EOS hydro

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

## Phase 5 — one-dimensional reactive flow

- [x] conserved `rho*Y_k` state coupled to NASA7 thermodynamics
- [x] temperature recovery from conserved total energy
- [x] general-EOS physical and Rusanov fluxes
- [x] exact species-flux closure to total mass flux
- [x] PCM lower-order baseline
- [x] frozen-composition characteristic PLM with MC/minmod limiting
- [x] composition-dependent CFL timestep
- [x] periodic and outflow boundaries
- [x] cell-local adiabatic constant-volume chemistry
- [x] reaction-hydro-reaction Strang splitting
- [x] homogeneous-field reduction to the 0D reactor
- [x] smooth entropy-wave second-order convergence
- [x] nonuniform reactive-hotspot structural regression
- [x] high-resolution hotspot comparison against PCM
- [ ] pressure-dependent/full H2/O2 chemistry in CFD
- [ ] molecular transport
- [ ] multidimensional reactive flow

## Verified `0.10.0` local results

| Gate | Result |
|---|---:|
| Complete Debug suite | `55/55` passed |
| Complete Release suite | `55/55` passed |
| Reactive entropy-wave density L1, 40 cells | `1.51594309e-4` |
| Reactive entropy-wave density L1, 80 cells | `3.43297270e-5` |
| Reactive entropy-wave density L1, 160 cells | `7.01334896e-6` |
| Observed order, 40 to 80 | `2.142685` |
| Observed order, 80 to 160 | `2.291283` |
| Hotspot completed steps | `475` |
| Hotspot maximum conservation error | `2.03e-16` |
| Hotspot pressure span | `1.0382765e3 Pa` |
| Hotspot maximum velocity magnitude | `2.2117421e1 m/s` |
| Hotspot temperature interval | `1328.97 to 1545.46 K` |
| Hotspot PLM normalized L1, 32 cells | `6.22861249e-2` |
| Hotspot PLM normalized L1, 64 cells | `1.46900072e-2` |
| Hotspot PCM normalized L1, 32 cells | `2.02487281e-1` |
| Hotspot PCM normalized L1, 64 cells | `1.54357039e-1` |

The existing `0.9.0` Cantera trajectory and exact-state production-rate comparisons remain active when `PELEF_ENABLE_CANTERA_REFERENCE=ON`.

## Next implementation slice

Move the pressure-dependent complete H2/O2 mechanism and an implicit cell reactor onto the verified Strang-split flow interface. Keep the four-reaction path as a regression baseline. After that, add molecular viscosity, thermal conduction, and species diffusion before extending reactive flow to two dimensions.
