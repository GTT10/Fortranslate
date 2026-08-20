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
- [x] periodic two-dimensional general-EOS CTU hydro
- [x] directional general-EOS HLLC with momentum rotation
- [x] full-state transverse correction with EOS positivity limiting

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
- [x] general-EOS physical, Rusanov, and HLLC fluxes
- [x] exact species-flux closure to total mass flux
- [x] PCM lower-order baseline
- [x] frozen-composition characteristic PLM with MC/minmod limiting
- [x] monotone primitive PPM with fourth-order face interpolation and SSPRK3
- [x] PeleC-style five-point characteristic PPM profile integration
- [x] species/transverse-velocity tracing on the middle wave
- [x] optional PeleC one-dimensional shock flattening
- [x] optional bounded Colella--Woodward contact steepening
- [x] composition-dependent CFL timestep
- [x] periodic and outflow boundaries
- [x] cell-local adiabatic constant-volume chemistry
- [x] reaction-hydro-reaction Strang splitting
- [x] homogeneous-field reduction to the 0D reactor
- [x] smooth entropy-wave second-order convergence
- [x] smooth composition-wave second-order convergence
- [x] discontinuous material-contact HLLC/Rusanov comparison
- [x] PPM/PLM material-contact resolution comparison
- [x] smooth PPM entropy/composition-wave convergence
- [x] smooth characteristic-PPM entropy/composition-wave convergence
- [x] characteristic-PPM/contact-steepening material-interface comparison
- [x] periodic reactive-mixture shock-flattening regression
- [x] nonuniform reactive-hotspot structural regression
- [x] nonuniform reactive-hotspot HLLC application regression
- [x] high-resolution hotspot comparison against PCM
- [ ] pressure-dependent/full H2/O2 chemistry in CFD
- [ ] molecular transport
- [x] periodic two-dimensional reactive CTU flow

## Phase 6 — two-dimensional reactive flow

- [x] separate `pelef_reactive_2d` driver and namelist
- [x] periodic uniform Cartesian reactive state and temperature field
- [x] x/y directional Rusanov and HLLC fluxes
- [x] frozen-composition characteristic PLM in both directions
- [x] conservative provisional/final CTU face fluxes
- [x] full-state transverse correction with bisection positivity limiting
- [x] species closure preserved through the same correction as mass and energy
- [x] composition-dependent two-dimensional CFL timestep
- [x] cell-local chemistry and Strang splitting
- [x] exact oblique entropy-wave convergence
- [x] transverse-correction signature and no-degradation gate
- [x] one-dimensional dimensional-reduction gate
- [x] periodic reactive-mixture vortex regression
- [x] two-dimensional reacting-hotspot application regression
- [ ] physical wall, inflow, and outflow boundaries
- [ ] multidimensional characteristic PPM and transverse PPM tracing
- [ ] molecular transport

## Verified `0.14.0` results

| Gate | Result |
|---|---:|
| Complete local Release suite without optional Cantera | `74/74` passed |
| New reactive-2D Debug tests | `6/6` passed |
| GitHub Actions Debug with Cantera | `75/75` target gate |
| GitHub Actions Release with Cantera | `75/75` target gate |
| Diagonal-wave density L1, 12 x 12 | `2.09551654e-4` |
| Diagonal-wave density L1, 24 x 24 | `7.15524055e-5` |
| Diagonal-wave density L1, 48 x 48 | `1.61190312e-5` |
| Diagonal-wave observed orders | `1.550234`, `2.150235` |
| 24 x 24 density L1 without CTU correction | `7.37180810e-5` |
| 1D/2D dimensional-reduction relative difference | below `3e-12` |
| Vortex maximum conservation error, 48 x 48 | `1.97e-15` |
| Vortex pressure span, 48 x 48 | `3.619e2 Pa` |
| 2D hotspot completed steps | `59` |
| 2D hotspot maximum conservation error | `9.31e-16` |
| 2D hotspot pressure span | `2.88685 Pa` |
| 2D hotspot maximum velocity | `6.75544e-3 m/s` |
| 2D hotspot temperature span | `2.44044e2 K` |
| 2D hotspot maximum H2O mass fraction | `1.53117e-4` |
| 2D hotspot maximum OH mass fraction | `1.58377e-5` |

## Verified `0.13.0` results

| Gate | Result |
|---|---:|
| Complete local Debug suite without optional Cantera | `68/68` passed |
| Complete local Release suite without optional Cantera | `68/68` passed |
| GitHub Actions Debug with Cantera | `69/69` passed |
| GitHub Actions Release with Cantera | `69/69` passed |
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
| Hotspot PPM normalized L1, 32 cells | `3.35268374e-2` |
| Hotspot PPM normalized L1, 64 cells | `5.08432812e-3` |
| HLLC composition-wave H2 L1, 40 cells | `2.73882775e-5` |
| HLLC composition-wave H2 L1, 80 cells | `6.85214336e-6` |
| HLLC composition-wave H2 L1, 160 cells | `1.65453783e-6` |
| HLLC composition-wave observed order, 40 to 80 | `1.998931` |
| HLLC composition-wave observed order, 80 to 160 | `2.050127` |
| PPM entropy-wave density L1, 32 cells | `3.02675468e-4` |
| PPM entropy-wave density L1, 64 cells | `7.43113334e-5` |
| PPM entropy-wave density L1, 128 cells | `1.75776010e-5` |
| PPM entropy-wave observed orders | `2.026118`, `2.079844` |
| PPM composition-wave H2 L1, 32 cells | `3.14258101e-5` |
| PPM composition-wave H2 L1, 64 cells | `7.19623593e-6` |
| PPM composition-wave H2 L1, 128 cells | `1.61225051e-6` |
| PPM composition-wave observed orders | `2.126636`, `2.158167` |
| Characteristic-PPM entropy-wave density L1, 32 cells | `2.49716040e-4` |
| Characteristic-PPM entropy-wave density L1, 64 cells | `4.88601959e-5` |
| Characteristic-PPM entropy-wave density L1, 128 cells | `9.87993492e-6` |
| Characteristic-PPM entropy-wave observed orders | `2.353557`, `2.306086` |
| Characteristic-PPM composition-wave H2 L1, 32 cells | `3.31512363e-5` |
| Characteristic-PPM composition-wave H2 L1, 64 cells | `7.50581812e-6` |
| Characteristic-PPM composition-wave H2 L1, 128 cells | `1.45326368e-6` |
| Characteristic-PPM composition-wave observed orders | `2.142981`, `2.368713` |
| Material-contact HLLC H2 L1 | `1.14269289e-4` |
| Material-contact PPM/HLLC H2 L1 | `8.75871389e-5` |
| Material-contact characteristic-PPM/HLLC H2 L1 | `7.35878653e-5` |
| Material-contact bounded-steepening H2 L1 | `2.50998077e-5` |
| Material-contact Rusanov H2 L1 | `1.87834655e-4` |
| Material-contact HLLC improvement factor | `1.6438` |
| Periodic shock pressure overshoot | `0` |
| Periodic shock flattening state-difference signature | `8.08412674e1` |
| Hotspot characteristic-PPM normalized L1, 32 cells | `6.00777009e-2` |
| Hotspot characteristic-PPM normalized L1, 64 cells | `3.07014574e-2` |

The existing `0.9.0` Cantera trajectory and exact-state production-rate comparisons remain active when `PELEF_ENABLE_CANTERA_REFERENCE=ON`.

## Next implementation slice

Add multidimensional characteristic PPM normal prediction and transverse PPM corrections to the reactive 2D path. In parallel, introduce mixture viscosity, thermal conduction, and species diffusion before moving the pressure-dependent complete H2/O2 mechanism and an implicit cell reactor onto the CFD interface.
