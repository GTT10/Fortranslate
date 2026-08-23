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
- [x] time-traced characteristic PPM normal prediction in x and y

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
- [x] third-body efficiencies
- [x] pressure falloff and Troe forms
- [ ] direct Cantera/CHEMKIN YAML parser
- [x] generated concentration and mass-fraction Jacobians
- [x] adaptive implicit backward Euler with step-doubling error control
- [x] complete ten-species, 29-reaction H2/O2 mechanism
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
- [x] pressure-dependent/full H2/O2 chemistry in CFD
- [x] qualified dilute-gas molecular transport in 1D
- [x] viscous, conductive, barodiffusive, correction-velocity, and enthalpy fluxes
- [x] explicit SSPRK2 transport and parabolic timestep gate
- [x] periodic transport-pulse and analytical shear-wave regressions
- [x] periodic two-dimensional reactive CTU flow

## Phase 6 — two-dimensional reactive flow

- [x] separate `pelef_reactive_2d` driver and namelist
- [x] periodic uniform Cartesian reactive state and temperature field
- [x] x/y directional Rusanov and HLLC fluxes
- [x] frozen-composition characteristic PLM in both directions
- [x] time-traced frozen-composition characteristic PPM in both directions
- [x] optional PeleC shock flattening on both normal predictors
- [x] optional bounded contact steepening on both normal predictors
- [x] conservative provisional/final CTU face fluxes
- [x] full-state transverse correction with bisection positivity limiting
- [x] species closure preserved through the same correction as mass and energy
- [x] composition-dependent two-dimensional CFL timestep
- [x] cell-local chemistry and Strang splitting
- [x] exact oblique entropy-wave convergence
- [x] exact oblique composition-wave convergence
- [x] transverse-correction signature and no-degradation gate
- [x] x-normal and y-normal one-dimensional dimensional-reduction gates
- [x] characteristic-PPM cell-level flattening/steepening signatures
- [x] oblique pressure-ratio-three shock positivity/flattening regression
- [x] two-dimensional material-contact sharpening regression
- [x] periodic reactive-mixture vortex regression
- [x] two-dimensional reacting-hotspot application regression
- [x] characteristic-PPM two-dimensional reacting-hotspot regression
- [x] slip/no-slip physical walls
- [x] adiabatic/isothermal wall temperature conditions
- [x] fixed-state inflow and zero-gradient outflow
- [ ] complete PeleC multidimensional PPM transverse/corner tracing
- [x] molecular transport


## Phase 7 — molecular transport

- [x] seven-species Lennard-Jones transport database pinned to Cantera
- [x] Chapman--Enskog pure viscosity and binary diffusion
- [x] Wilke mixture viscosity
- [x] modified-Eucken pure conductivity and Mathur mixture conductivity
- [x] mixture-averaged species diffusion
- [x] one-dimensional Newtonian viscous stresses and viscous energy work
- [x] Fourier heat conduction
- [x] mole-fraction diffusion with optional barodiffusion
- [x] correction velocity enforcing zero net diffusive mass flux
- [x] species-enthalpy diffusion contribution to total-energy flux
- [x] explicit SSPRK2 transport advance and parabolic CFL limit
- [x] symmetric chemistry/transport/hydro operator composition
- [x] Cantera 3.2 coefficient qualification probe
- [x] second-order viscous shear-wave convergence
- [x] species and temperature smoothing with periodic conservation
- [ ] Soret and Dufour effects
- [ ] multicomponent Stefan--Maxwell diffusion
- [ ] full PelePhysics polynomial/polar transport parity
- [x] molecular transport in the reactive 2D CTU path

## Verified `0.16.0` results

| Gate | Result |
|---|---:|
| Complete local Debug suite without optional Cantera | `86/86` passed |
| Complete local Release suite without optional Cantera | `86/86` passed |
| GitHub Actions Debug with Cantera | `87/87` target gate |
| GitHub Actions Release with Cantera | `87/87` target gate |
| Viscous shear-wave density grids | `32`, `64`, `128` |
| Viscous shear-wave transverse-velocity L1 | `2.352743e-6`, `5.849738e-7`, `1.460204e-7` |
| Viscous shear-wave observed orders | `2.007900`, `2.002202` |
| 1000 K stoichiometric mixture viscosity | `4.19833898e-5 Pa s` |
| 1000 K stoichiometric mixture conductivity | `1.24288597e-1 W/(m K)` |
| 1000 K H2/O2/N2 mixture diffusion coefficients | `8.12964482e-4`, `1.98336193e-4`, `1.80234445e-4 m2/s` |
| 96-cell transport-pulse steps | `297` |
| Transport-pulse H2 range reduction | `5.27325851e-3` to `5.10182280e-3` |
| Transport-pulse maximum species-mass error | `1.57e-15` |
| Transport-pulse maximum conservation error | `1.73e-18` |
| Transport-pulse pressure span | `29.0691 Pa` |
| Transport-pulse maximum velocity | `6.28770e-2 m/s` |

## Verified `0.15.0` results

| Gate | Result |
|---|---:|
| Complete local Debug suite without optional Cantera | `80/80` passed |
| Complete local Release suite without optional Cantera | `80/80` passed |
| GitHub Actions Debug with Cantera | `81/81` target gate |
| GitHub Actions Release with Cantera | `81/81` target gate |
| Characteristic-PPM diagonal density L1, 16 x 16 | `1.44130250e-4` |
| Characteristic-PPM diagonal density L1, 32 x 32 | `3.80399208e-5` |
| Characteristic-PPM diagonal density L1, 64 x 64 | `9.10858005e-6` |
| Characteristic-PPM density observed orders | `1.921787`, `2.062216` |
| 32 x 32 density L1 without CTU correction | `4.33952478e-5` |
| Characteristic-PPM H2 L1, 16 x 16 | `3.74405262e-5` |
| Characteristic-PPM H2 L1, 32 x 32 | `9.88530404e-6` |
| Characteristic-PPM H2 L1, 64 x 64 | `2.36662586e-6` |
| Characteristic-PPM H2 observed orders | `1.921243`, `2.062454` |
| x/y 1D dimensional-reduction relative difference | below `3e-12` |
| 2D material-contact H2 L1 without steepening | `2.65683522e-4` |
| 2D material-contact H2 L1 with bounded steepening | `1.66701302e-4` |
| Oblique 2D shock pressure overshoot | `0` |
| Oblique 2D shock flattening state-difference signature | `1.55397791e2` |
| Characteristic-PPM 2D hotspot completed steps | `59` |
| Characteristic-PPM 2D hotspot maximum conservation error | `2.40e-15` |
| Characteristic-PPM 2D hotspot pressure span | `2.82161 Pa` |
| Characteristic-PPM 2D hotspot maximum velocity | `6.91813e-3 m/s` |
| Characteristic-PPM 2D hotspot temperature span | `2.44044e2 K` |
| Characteristic-PPM 2D hotspot maximum H2O mass fraction | `1.53117e-4` |
| Characteristic-PPM 2D hotspot maximum OH mass fraction | `1.58377e-5` |

## Verified `0.14.0` results

| Gate | Result |
|---|---:|
| Complete local Release suite without optional Cantera | `74/74` passed |
| Complete local Debug suite without optional Cantera | `74/74` passed |
| GitHub Actions Debug with Cantera | `75/75` passed |
| GitHub Actions Release with Cantera | `75/75` passed |
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

Add catalytic or prescribed-species wall fluxes and characteristic open-boundary
handling after the present slip/no-slip, thermal-wall, fixed-inflow, and extrapolated-
outflow milestone. Complete PeleC multidimensional PPM corner coupling remains a
separate parity task.


## Verified `0.17.0` results

| Gate | Result |
|---|---:|
| Local Debug suite | `90/90` passed |
| Local Release suite | `90/90` passed |
| 2D shear-wave L1, 16 x 16 | `2.603378e-6` |
| 2D shear-wave L1, 32 x 32 | `6.316041e-7` |
| 2D shear-wave L1, 64 x 64 | `1.563289e-7` |
| Observed shear-wave orders | `2.043293`, `2.014436` |
| Reacting transport-hotspot steps | `49` |
| Reacting transport-hotspot conservation error | `1.33063574e-15` |
| Reacting transport-hotspot pressure span | `8.34842515 Pa` |
| Reacting transport-hotspot maximum speed | `9.62290759e-2 m/s` |
| Reacting transport-hotspot temperature span | `2.40627606e2 K` |


## Physical boundaries (`0.18.0`)

- [x] typed x/y boundary-face configuration
- [x] matched periodic-pair validation
- [x] moving slip/no-slip walls
- [x] adiabatic/isothermal wall heat flux
- [x] zero solid-wall species diffusion
- [x] fixed inflow and extrapolated outflow
- [x] boundary-aware PLM/PPM sampling and explicit face divergence
- [x] Couette, thermal-wall, and inflow/outflow regressions
- [ ] catalytic walls, NSCBC, and characteristic non-reflecting outflow

## Verified `0.18.0` results

| Gate | Result |
|---|---:|
| Complete local Release suite | `92/92` passed |
| Complete local Debug suite | `92/92` passed |
| Couette upper-wall velocity | `20 m/s` |
| Couette maximum cell-center velocity | `19.6875 m/s` |
| Uniform inflow speed | `75 m/s` |


## Pressure-dependent chemistry in CFD (`0.19.0`)

- [x] ten-species thermodynamics and transport
- [x] 29 reversible reactions with third-body, falloff, and Troe forms
- [x] implicit constant-volume chemistry in 1D and 2D cells
- [x] runtime elementary/full mechanism selection
- [x] uniform 0D/1D/2D reduction gate


## MPI domain decomposition (`0.20.0`–`0.24.0`)

- [x] uneven contiguous 1D decomposition
- [x] nonblocking periodic halo exchange
- [x] global CFL and conservation reductions
- [x] ordered `MPI_Gatherv` output
- [x] 15-component halo verification
- [x] 1/2/4/8-rank field parity
- [x] distributed conservative multispecies hydro
- [x] distributed general-EOS reactive hydro
- [x] distributed molecular transport
- [x] distributed implicit chemistry scheduling
- [x] globally synchronized adaptive step rejection and rollback
- [x] coupled reaction--transport--hydro--transport--reaction splitting

The `0.24.0` milestone completes the planned one-dimensional MPI slice. AMR,
multidimensional MPI decomposition, load balancing, and accelerator execution
remain later porting-plan phases rather than claims of this milestone.

## AMR foundation (`0.25.0`)

- [x] typed one-dimensional level and patch/box metadata
- [x] static two-level hierarchy with a strictly nested fine patch
- [x] integer refinement ratios and fine/coarse spacing relationship
- [x] MC-limited conservative piecewise-linear prolongation
- [x] volume-average restriction and average-down synchronization
- [x] refinement-ratio level subcycling schedule
- [x] time-integrated coarse/fine flux register
- [x] conservative two-interface reflux
- [x] composite-integral conservation gate
- [x] solution-driven tagging (`0.26.0`)
- [x] dynamic regridding (`0.26.0`)
- [ ] arbitrary multiple levels
- [x] AMR reactive-flow application (`0.27.0`)
- [x] AMR molecular transport (`0.29.0`)

The foundation is state-width independent and intentionally serial. MPI patch
distribution follows only after the AMR ownership and regridding model is
qualified.

## Dynamic AMR regridding (`0.26.0`)

- [x] component-selectable relative-gradient tagging
- [x] absolute gradient floor and local normalization floor
- [x] deterministic tag buffering and minimum patch width
- [x] explicit rejection of unsupported physical-boundary tags
- [x] conservative fine-patch creation, relocation, resizing, and removal
- [x] average-down before an old refined region is discarded
- [x] exact transfer of overlapping fine cells at unchanged refinement ratio
- [x] composite-integral invariance across regrid operations

The current planner deliberately produces one strictly interior patch. It does
not yet cluster disjoint tag sets, refine physical boundaries, or construct
more than two levels.

## Reactive AMR application (`0.27.0`)

- [x] executable two-level reactive 1D driver and namelist configuration
- [x] general-EOS coarse/fine conserved states and temperatures
- [x] refinement-ratio fine hydro subcycling
- [x] time-interpolated coarse data for fine-patch ghost cells
- [x] coarse/fine interface flux accumulation and reflux every coarse step
- [x] covered-cell average-down after hydro and chemistry
- [x] hierarchy-wide chemistry--hydro--chemistry Strang composition
- [x] transactional rollback of coarse and fine states on failed updates
- [x] periodic solution-driven patch creation, movement, resizing, or removal
- [x] composite CSV output with ordered exact domain coverage
- [x] reacting-hotspot positivity, closure, synchronization, and conservation

This milestone initially qualified first-order PCM hydro with elementary or
full-H2/O2 chemistry.

## Reactive AMR PLM (`0.28.0`)

- [x] MC/minmod-limited primitive-variable coarse and fine reconstruction
- [x] face-state density and pressure positivity fallback
- [x] nonnegative species-face normalization and closure
- [x] SSPRK2 level advancement
- [x] time-averaged SSPRK2 interface fluxes supplied to reflux
- [x] periodic physical-boundary flux identity
- [x] midpoint coarse-time ghosts for each fine PLM substep
- [x] moving-contact error lower than the matching PCM AMR run
- [x] composite conservation retained with PLM and dynamic regrid

Characteristic PPM/WENO AMR reconstruction, multiple patches, more than two
levels, and MPI patch ownership remain open.

## Reactive AMR molecular transport (`0.29.0`)

- [x] viscosity, Fourier conduction, and mixture-averaged species diffusion
- [x] optional barodiffusion, correction velocity, and species-enthalpy flux
- [x] parabolic coarse/fine timestep limits
- [x] `r^2` fine transport subcycling with coarse-time ghost interpolation
- [x] coarse/fine center-distance gradients at patch interfaces
- [x] SSPRK2 stage-averaged diffusive flux registration
- [x] diffusive reflux and covered-cell average-down after each half step
- [x] reaction--transport--hydro--transport--reaction AMR composition
- [x] transactional rollback across every split operator
- [x] periodic all-component and per-species conservation gates
- [x] conduction/reference smoothing and final synchronization gates

Soret, Dufour, multicomponent Stefan--Maxwell diffusion, PelePhysics transport
polynomial parity, multiple patches, and arbitrary multilevel recursion remain
outside this milestone.
