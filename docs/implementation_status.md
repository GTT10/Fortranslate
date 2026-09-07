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
- [x] three-dimensional periodic PCM/SSPRK2 Euler application
- [x] x/y/z whole-step reduction and smooth-wave conservation/convergence

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
- [x] qualified general-EOS PeleC-style Riemann parity
- [x] periodic two-dimensional general-EOS CTU hydro
- [x] directional general-EOS HLLC with momentum rotation
- [x] periodic three-dimensional general-EOS multispecies PCM/PLM SSPRK2 hydro
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
- [x] pinned build-time Cantera YAML ingestion for the qualified H2/O2 subset
- [ ] general Cantera YAML and CHEMKIN mechanism ingestion
- [x] generated concentration and mass-fraction Jacobians
- [x] adaptive implicit backward Euler with step-doubling error control
- [x] complete ten-species, 29-reaction H2/O2 mechanism
- [ ] detailed hydrocarbon mechanism

## Phase 5 — one-dimensional reactive flow

- [x] conserved `rho*Y_k` state coupled to NASA7 thermodynamics
- [x] temperature recovery from conserved total energy
- [x] general-EOS physical, Rusanov, HLLC, and PeleC-style fluxes
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
- [x] x/y directional Rusanov, HLLC, and PeleC-style fluxes
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
- [x] prescribed zero-net-mass wall species and enthalpy fluxes
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

Add catalytic surface-rate models on top of the prescribed wall-flux contract,
then characteristic open-boundary handling. Complete PeleC multidimensional
PPM corner coupling remains a separate parity task.


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
- [x] prescribed zero-net-mass wall species and enthalpy flux
- [ ] catalytic surface kinetics, NSCBC, and characteristic non-reflecting outflow

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
- [x] arbitrary multiple levels in hierarchy primitives (`0.30.0`)
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

Characteristic PPM/WENO AMR reconstruction, multiple patches, dynamic
regridding beyond two levels, and MPI patch ownership remain open.

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

## Arbitrary-depth AMR hierarchy foundation (`0.30.0`)

- [x] allocatable sequence of adjacent coarse/fine level interfaces
- [x] contiguous level numbering and strict physical nesting validation
- [x] independent integer refinement ratio at every interface
- [x] cumulative cell spacing, physical bounds, and subcycle products
- [x] variable-size state storage for every level
- [x] recursive conservative prolongation from the root level
- [x] deepest-to-root restriction and covered-cell average-down
- [x] one flux register per adjacent level pair
- [x] deepest-to-root reflux followed by average-down
- [x] arbitrary-depth composite integral over uncovered cells
- [x] four-level mixed-ratio conservation and synchronization gate

The hierarchy primitives have no fixed level-count limit. The runnable reactive
application still owns two levels; the static recursive reactive engine below
uses this hierarchy while recursive regridding and output remain separate.

## Arbitrary-depth reactive AMR advancement (`0.31.0`)

- [x] allocatable conserved-state and temperature fields at every level
- [x] root physical ghosts and coarse/fine ghosts at every nested interface
- [x] global CFL limit reduced from every level into a root time step
- [x] recursive `r`-way hydro subcycling with coarse-time interpolation
- [x] recursive `r^2`-way molecular-transport subcycling
- [x] one time-integrated flux register for every active parent/child advance
- [x] deepest-to-root reflux, average-down, and temperature recovery
- [x] hierarchy-wide chemistry and `R-T-H-T-R` operator composition
- [x] whole-hierarchy transactional rollback on a failed operator
- [x] arbitrary-depth composite conserved integrals
- [x] three-level transport/chemistry conservation and positivity regression

The recursive engine accepts runtime patch bounds and refinement ratios with no
fixed depth. Dynamic tag-driven construction, regridding during time advance,
composite CSV output, multiple patches, physical-boundary refinement, and MPI
patch distribution are not claimed by this milestone.

## Dynamic multilevel reactive AMR application (`0.32.0`)

- [x] `amr_max_levels` runtime configuration with a two-level-compatible default
- [x] recursive normalized-gradient tag planning at every active depth
- [x] strict nesting through suppression of interior patch-edge tags
- [x] automatic termination when a parent has no refinement tags
- [x] deepest-to-root average-down before any hierarchy replacement
- [x] conservative prolongation of every rebuilt child patch
- [x] unchanged-plan detection that retains the complete existing hierarchy
- [x] periodic multilevel regrid evaluation in the simulation loop
- [x] exact final-time clipping and five-component composite diagnostics
- [x] spatially ordered recursive composite CSV output
- [x] executable three-level hotspot case and output coverage checker
- [x] forced hierarchy-change and regrid-conservation regression

When a multilevel plan changes, this milestone deliberately transfers the old
composite through the synchronized root and reconstructs children. This
preserves conserved integrals but not old fine-scale values in overlap regions.
Multiple patches, physical-boundary refinement, MPI patch ownership, and
load-balanced distribution remain open.

## Multilevel overlap-preserving regrid transfer (`0.33.0`)

- [x] old/new physical-bound intersection at every common fine level
- [x] spacing and cell-boundary alignment checks before direct transfer
- [x] exact conserved-state and temperature copy for aligned overlap cells
- [x] safe conservative-prolongation fallback when level spacing changes
- [x] deepest-to-root synchronization after all overlap copies
- [x] cumulative transferred-cell diagnostics in the executable
- [x] forced hierarchy-change regression with bitwise deepest-level retention
- [x] simultaneous composite-conservation gate across the changed hierarchy

The qualified transfer assumes one Cartesian patch per level. It maps aligned
cells at equal level spacing; arbitrary refinement-ratio remapping, multiple
patch intersection, physical-boundary refinement, and distributed ownership
remain separate.

## Multilevel reactive AMR characteristic PPM (`0.34.0`)

- [x] four exterior conserved-state and temperature layers per level
- [x] periodic/outflow physical PPM ghost construction
- [x] conservative MC-limited parent-to-fine subcell interpolation
- [x] parent start/end time interpolation at every fine hydro subcycle
- [x] primitive PPM and time-traced frozen-composition characteristic PPM
- [x] optional bounded contact steepening and shock flattening
- [x] SSPRK3 advancement on every recursively subcycled level
- [x] SSPRK3 effective face flux supplied to every flux register
- [x] PPM-aware strict-nesting buffer in dynamic tag planning
- [x] three-level conservation, synchronization, positivity, and closure gate
- [x] runnable dynamic three-level characteristic-PPM hotspot case

The qualified path holds the midpoint parent ghost interpolation fixed through
the three SSPRK stages of one fine substep. It retains one interior parent cell
beyond the four-layer footprint for limited spatial interpolation.

## Multilevel reactive AMR hybrid WENO5 (`0.35.0`)

- [x] PeleC-style five-point WENO5-JS edge reconstruction
- [x] PeleC-style five-point WENO5-Z edge reconstruction
- [x] exact constant/linear reproduction and fixed formula-parity points
- [x] optional hybrid replacement inside characteristic PPM
- [x] four-layer physical and coarse/fine ghost reuse at every level
- [x] recursive SSPRK3 advancement and effective-flux reflux reuse
- [x] three-level JS/Z conservation, synchronization, positivity, and closure
- [x] runnable dynamic three-level WENO5-Z hotspot case

This is a one-dimensional multilevel reactive qualification. WENO3-Z,
WENO7-Z, regular-grid two-dimensional WENO, physical-boundary refinement,
multiple patches, and distributed patch ownership remain separate.

## Complete PeleC hybrid-WENO scheme set (`0.36.0`)

- [x] WENO7-Z seven-point reconstruction (`amr_weno_scheme = 2`)
- [x] WENO3-Z three-point reconstruction (`amr_weno_scheme = 3`)
- [x] exact constant/linear reproduction and fixed formula-parity points
- [x] four-scheme fixed three-level conservation and synchronization gate
- [x] public dynamic three-level WENO7-Z and WENO3-Z hotspot cases

Together with the 0.35.0 WENO5-JS/Z slice, the four hybrid schemes currently
selected by PeleC `PPM.cpp` are available. Regular-grid two-dimensional WENO,
physical-boundary refinement, multiple patches, and distributed patch
ownership remain separate.

## Outflow physical-boundary refinement (`0.37.0`)

- [x] hierarchy geometry accepts left-, right-, or both-side parent contact
- [x] conservative endpoint prolongation without parent out-of-bounds access
- [x] one-sided reflux only at an existing coarse/fine interface
- [x] composite integration with zero uncovered cells on either side
- [x] fine face and four-layer PPM/WENO outflow ghost construction
- [x] recursive nested boundary contact at two successive interfaces
- [x] WENO7-Z conservation, synchronization, positivity, and closure gate
- [x] runnable tag-driven three-level boundary-hotspot case

The qualified dynamic path is outflow-only. A periodic child may wrap itself
only when it covers its complete parent. A single patch crossing the periodic
seam, multiple patches, and distributed patch ownership remain separate.

## Two-level multipatch AMR foundation (`0.38.0`)

- [x] ordered non-overlapping patch-set geometry over one parent level
- [x] boundary-contact and empty-set representations
- [x] set-wide conservative prolongation and average-down
- [x] one flux register per patch and transactional set-wide reflux
- [x] composite integration excluding every covered parent interval once
- [x] disconnected-tag clustering and adjacent-candidate coalescing
- [x] conservative patch-set creation, movement, repartition, and removal
- [x] exact same-resolution fine overlap transfer across old/new patch pairs
- [x] fixed two-level reactive WENO7-Z subcycling on two fine patches
- [x] conservation, synchronization, positivity, and species-closure gate

The `0.38.0` reactive qualification is hydro-only and uses separated patches,
so each coarse/fine side is a genuine interface and no same-level ghost
exchange is required.

## Two-level multipatch reactive physics (`0.39.0`)

- [x] hierarchy-wide chemistry half steps with set-wide average-down
- [x] per-patch SSPRK2 molecular-transport flux registers
- [x] `r^2` fine transport subcycling with parent-time ghost interpolation
- [x] combined hyperbolic/parabolic root timestep reduction across all patches
- [x] transactional `R-T-H-T-R` composition and missing-database rejection
- [x] reacting-hotspot conservation, synchronization, positivity, and closure

The fixed two-level separated-patch path now supports chemistry and molecular
transport. At the `0.39.0` milestone, the primary tag-driven application still
owned one patch per level; Decision 0032 subsequently integrates the two-level
patch-set engine into that application.

## Dynamic two-level multipatch application (`0.40.0`)

- [x] runtime multipatch selection and maximum tag-gap configuration
- [x] disconnected-tag clustering into an empty or separated patch collection
- [x] periodic-seam-safe and outflow-boundary-safe patch planning
- [x] interval-driven patch creation, movement, repartition, and removal
- [x] exact transfer of every aligned old/new fine-cell intersection
- [x] regrid-evaluation, accepted-regrid, and overlap-transfer counters
- [x] ordered exact-coverage multipatch composite CSV output
- [x] public chemistry-and-transport entropy-wave application gate
- [x] deterministic empty-to-two-to-moved-to-empty lifecycle regression

The public reactive AMR executable selects the dynamic two-level patch-set
engine when `amr_multipatch_enabled` is true. Static arbitrary-depth
multipatch recursion is qualified separately in `0.43.0`. Dynamic patch-tree
rebuilds, independently owned adjacent boxes and same-level exchange,
periodic-seam splitting, load balancing, and distributed patch ownership
remain separate.

## Arbitrary-depth multipatch tree foundation (`0.41.0`)

- [x] explicit parent ownership for every child patch
- [x] zero-or-more separated child sets per parent patch
- [x] deterministic parent-local to flattened-level child indexing
- [x] arbitrary relation count with mixed per-level refinement ratios
- [x] physical-extent and parent-ownership hierarchy validation
- [x] root-to-leaf conservative patch-set prolongation
- [x] transactional deepest-to-root patch-set average-down
- [x] exact branching-tree composite integration
- [x] four-level `1/2/3/2`-patch mixed-ratio unit gate

This milestone qualifies static geometry and conservative field operations.
The existing arbitrary-depth reactive time integrator still owns one patch per
level. Patch-tree subcycling, reflux, chemistry/transport integration, dynamic
tagging and overlap transfer, adjacent boxes, same-level exchange, and
distributed ownership remain separate.

## Arbitrary-depth multipatch tree synchronization (`0.42.0`)

- [x] one flux-register array per parent-owned child set
- [x] exact nested register shape validation, including empty child sets
- [x] transactional field and register rollback across the full tree
- [x] deepest-to-root per-parent reflux followed by average-down
- [x] register reset after accepted synchronization
- [x] branched deepest-level flux-mismatch conservation gate

This milestone supplied the synchronization primitive consumed by the
`0.43.0` recursive hydro path. Reactive transport still needs to accumulate
its time-integrated interface fluxes into the nested register layout. Dynamic
tree rebuilds, adjacent boxes, same-level exchange, load balancing, and
distributed ownership remain separate.

## Arbitrary-depth reactive patch-tree hydro (`0.43.0`)

- [x] reactive state, temperature, and wide-ghost ownership per tree patch
- [x] all-patch CFL reduction to a stable root-level interval
- [x] depth-first parent advance with refinement-ratio child subcycling
- [x] time-interpolated parent-to-child coarse/fine ghost fill
- [x] per-child coarse and fine time-integrated flux accumulation
- [x] per-parent reflux and covered-cell average-down after child subcycles
- [x] transactional rollback of the complete reactive tree
- [x] four-level branched PCM hydro conservation and synchronization gate
- [x] exact `1/4/12/16` level-advance accounting gate
- [x] temperature, pressure, species positivity, and closure gates

This milestone qualifies static, strictly interior, separated patch trees for
reactive PCM hydro. Chemistry and molecular transport composition, dynamic
tagging and overlap-preserving tree rebuilds, adjacent boxes and same-level
exchange, physical-boundary children, load balancing, and distributed patch
ownership remain separate.

## Reactive patch-tree chemistry splitting (`0.44.0`)

- [x] elementary chemistry advancement on every stored patch
- [x] identical physical reaction interval across all tree levels
- [x] deepest-to-root average-down after each reaction half-step
- [x] temperature recovery and ghost refresh after reaction synchronization
- [x] symmetric chemistry--hydro--chemistry public advance
- [x] whole-tree rollback across both reaction halves and recursive hydro
- [x] chemistry-enabled versus hydro-only species-evolution gate
- [x] composite mass, momentum, and total-energy conservation gate
- [x] post-reaction synchronization, positivity, and species-closure gates

This milestone composes elementary chemistry with the static patch-tree PCM
hydro path. Molecular transport, dynamic tagging and overlap-preserving tree
rebuilds, same-level exchange, physical-boundary children, load balancing, and
distributed patch ownership remain separate.

## Reactive patch-tree molecular transport (`0.45.0`)

- [x] all-patch hydro and parabolic timestep reduction
- [x] cumulative refinement-ratio-squared transport subcycling
- [x] midpoint time interpolation for parent-to-child transport ghosts
- [x] per-child coarse and fine diffusive flux accumulation
- [x] recursive transport reflux and covered-cell average-down
- [x] separate per-level transport-advance accounting
- [x] reaction--transport--hydro--transport--reaction public composition
- [x] four-level exact `2/16/96/256` transport-call gate
- [x] transport-enabled versus transport-disabled state-change gate
- [x] mass, momentum, energy, synchronization, and physical-state gates
- [x] missing transport database rejection without solution mutation

This milestone completes the physics composition on the static, strictly
interior patch-tree PCM path. Dynamic tagging and overlap-preserving tree
rebuilds, same-level exchange, physical-boundary children, load balancing, and
distributed patch ownership remain separate.

## Runtime reactive patch-tree rebuilds (`0.46.0`)

- [x] explicit runtime branching-plan replacement
- [x] identical-plan no-op detection and evaluation accounting
- [x] old-tree deepest-to-root synchronization before rebuild
- [x] conservative root-to-new-tree prolongation
- [x] physical-coordinate same-spacing overlap transfer at every common level
- [x] overlap transfer independent of changed parent ownership
- [x] final deepest-to-root synchronization and temperature recovery
- [x] time, hydro, transport, and regrid counter preservation
- [x] whole-tree rollback for invalid plans or failed transfer
- [x] moved four-level conservation and deepest exact-retention gate

This milestone qualifies plan-driven dynamic trees. Automatic tag clustering
at every parent is added in `0.47.0`, and same-level sibling exchange in
`0.48.0`; physical-boundary children, load balancing, and distributed
ownership remain separate.

## Tag-driven reactive patch-tree rebuilds (`0.47.0`)

- [x] synchronized-root canonical input for deterministic planning
- [x] normalized-gradient tagging on every prospective parent patch
- [x] parent-local disconnected-tag clustering and minimum-width expansion
- [x] stencil-support reservation for PCM/PLM and wide PPM ghosts
- [x] deterministic parent-order flattening into branching level plans
- [x] independent termination of untagged branches
- [x] recursive plan construction through `amr_max_levels`
- [x] tag-plan connection to the transactional overlap-preserving rebuild
- [x] root-only to four-level `1/2/2/2` branch-creation gate
- [x] unchanged-plan no-op and invalid-tag-request rollback gates

This milestone qualifies automatic interior patch-tree planning and rebuilds
for one common configured refinement ratio. Same-level adjacent-patch exchange
is added in `0.48.0`; physical-boundary children, load balancing, and
distributed patch ownership remain separate.

## Adjacent patch-tree same-level exchange (`0.48.0`)

- [x] patch-tree-only opt-in for independently owned adjacent child intervals
- [x] parent-local face ghost exchange by global fine-cell index
- [x] PPM/WENO four-layer exchange across adjacent sibling chains
- [x] one shared time-integrated flux per fine/fine face
- [x] conservative correction of both cells adjacent to the shared face
- [x] internal fine/fine register-side suppression before reflux
- [x] identical interface ownership for hydro and molecular transport
- [x] exact post-initialization and post-step ghost exchange gates
- [x] PPM hydro conservation and `1/4` subcycle-accounting gate
- [x] transport conservation and `2/16` parabolic-accounting gate
- [x] parent-child synchronization after adjacent hydro and transport

This milestone qualifies same-process sibling exchange for strictly interior
patch-tree children. Periodic-seam splitting, physical-boundary children,
load balancing, and distributed MPI patch ownership remain separate.

## MPI AMR patch distribution bridge (`0.49.0`)

- [x] collective validation of replicated patch-tree topology and extent
- [x] deterministic cell-weighted owner for every root/fine patch
- [x] exact one-owner accounting of patch and cell work
- [x] owner-authoritative full-patch field synchronization
- [x] parent-local adjacent-sibling discovery in global fine indices
- [x] owner-authoritative cross-rank halo exchange through four layers
- [x] identical collective ordering independent of local ownership
- [x] collective rejection of rank-inconsistent valid hierarchies
- [x] 1/2/4/8-rank Release and Debug gates

This milestone establishes MPI ownership metadata and the communication bridge
while keeping the serial AMR modules free of MPI. Patch field arrays remain
replicated on every rank. At `0.49.0`, owner-only physics, distributed
fine/fine flux reconciliation, sparse storage, regrid migration, and scalable
point-to-point exchange remained separate integration work.

## Owner-only MPI AMR chemistry (`0.50.0`)

- [x] full reactive state, temperature, and wide-ghost owner synchronization
- [x] collective reactive layout and state-width preflight
- [x] one chemistry integration on the exclusive owner of each patch
- [x] communicator-wide acceptance after every patch reactor call
- [x] owner broadcast after each accepted patch update
- [x] deepest-to-root average-down and temperature recovery
- [x] physical/coarse-fine/sibling ghost reconstruction after chemistry
- [x] exact all-rank rollback after a deep owner failure
- [x] exact global patch-call accounting independent of rank count
- [x] serial patch-tree reacting-field parity and conservation
- [x] 1/2/4/8-rank Release and Debug gates

The chemistry operator is the first patch-tree physics path executed only on
owners. Replicated storage remains intentional for this bridge. Owner-only
hydro and transport recursion, cross-rank shared-flux reconciliation, sparse
rank-local storage, and regrid migration remain the next distributed AMR work.

## Owner-only MPI AMR hydro (`0.51.0`)

- [x] shared serial/MPI one-patch PCM/PLM/PPM hydro kernel
- [x] identical collective parent/child/substep traversal on every rank
- [x] one finite-volume patch update on the exclusive owner per substep
- [x] collective patch acceptance and owner state/face-flux broadcasts
- [x] time-interpolated parent ghosts and adjacent-sibling halo replacement
- [x] cross-owner shared time-integrated fine/fine flux reconciliation
- [x] replicated coarse/fine flux registers, reflux, and average-down
- [x] exact global rollback after a deep owner hydro failure
- [x] exact owner and per-level cumulative subcycle accounting
- [x] serial parity and composite conservation for four-level PCM hydro
- [x] serial parity and conservation for six adjacent PPM children
- [x] 1/2/4/8-rank Release and Debug gates

The finite-volume update itself now executes only on patch owners, including
recursive level-ratio subcycling. The bridge broadcasts complete owner fluxes
and states, so synchronization and flux-register work still run on replicated
trees. Owner-only molecular transport, sparse rank-local allocation,
stage-synchronous point-to-point halos, and regrid migration remain pending.

## Owner-only MPI AMR molecular transport (`0.52.0`)

- [x] shared serial/MPI one-patch SSPRK2 molecular-transport kernel
- [x] owner-only viscous, conductive, and species-diffusion patch updates
- [x] exact cumulative `r^2` child subcycling on every relation
- [x] collective acceptance and owner state/effective-flux broadcasts
- [x] time-interpolated parent and adjacent-sibling transport ghosts
- [x] cross-owner shared time-integrated diffusive face fluxes
- [x] replicated diffusive registers, reflux, average-down, and recovery
- [x] exact global rollback after a deep owner transport failure
- [x] exact owner and per-level parabolic call accounting
- [x] serial parity and conservation for a four-level branched tree
- [x] serial parity and conservation for six adjacent children
- [x] 1/2/4/8-rank Release and Debug gates

All three patch-tree physics operators now have owner-only MPI entry points.
They remain separate qualified transactions over replicated trees. A combined
distributed full-physics transaction, sparse rank-local storage,
stage-synchronous point-to-point halos, and regrid migration remain pending.

## Transactional owner-only MPI AMR full physics (`0.53.0`)

- [x] distributed `R-T-H-T-R` composition using qualified owner operators
- [x] root-owner synchronization of time, step, counters, and regrid metadata
- [x] one identical outer backup on every rank before the first operator
- [x] communicator-wide acceptance after every component operator
- [x] zero accepted-call reporting for a rejected outer transaction
- [x] missing transport-database rejection before mutation
- [x] exact rollback after chemistry and transport precede hydro failure
- [x] exact global chemistry/hydro/transport call accounting
- [x] serial full-state, temperature, ghost, and bookkeeping parity
- [x] composite full-physics conservation on a four-level branched tree
- [x] 1/2/4/8-rank Release and Debug gates

The replicated MPI bridge now covers the complete patch-tree physics interval,
not only separate operators. Sparse rank-local patch allocation, ownership
migration after regrid, and scalable point-to-point stage communication remain
the next distributed-AMR architecture work.

## Sparse MPI AMR patch storage and owner migration (`0.54.0`)

- [x] replicated hierarchy and owner descriptors with rank-local payloads
- [x] all six reactive patch arrays absent on every non-owner rank
- [x] exact owner-authoritative replica-to-sparse scatter
- [x] exact sparse-to-replicated gather including global bookkeeping
- [x] global patch, cell, and field-value one-copy accounting
- [x] same-hierarchy owner-map migration with no retained non-owner payload
- [x] exact post-migration field, ghost, and bookkeeping reconstruction
- [x] 1/2/4/8-rank Release and Debug gates

This establishes the storage boundary needed for a truly distributed AMR
advance. The current physics entry points still use replicated trees, and
owner migration temporarily broadcasts one patch at a time. Direct sparse
physics, topology-changing regrid migration, and scalable point-to-point halo
and flux schedules remain pending.

## Direct chemistry on sparse MPI AMR patches (`0.55.0`)

- [x] chemistry integration only on locally allocated owner payloads
- [x] exact one-call-per-patch global accounting
- [x] deepest-to-root child streaming and owner-local average-down
- [x] owner-local temperature recovery after covered-parent replacement
- [x] root physical and parent/fine ghost refresh without a full replica
- [x] cross-owner adjacent normal and four-layer PPM ghost replacement
- [x] serial parity for four-level branched and adjacent-patch trees
- [x] exact sparse rollback after a deep owner chemistry failure
- [x] 1/2/4/8-rank Release and Debug gates

Chemistry is now the first AMR operator to consume sparse storage directly.
Synchronization currently broadcasts one child, parent, or sibling patch at a
time. Sparse hydro and transport, point-to-point schedules, and distributed
topology-changing regrid remain pending.

## Direct recursive hydro on sparse MPI AMR patches (`0.56.0`)

- [x] owner-only finite-volume updates on rank-local sparse payloads
- [x] mixed-ratio recursive level subcycling without a full tree replica
- [x] streamed parent start/end states for time-interpolated child ghosts
- [x] replicated compact coarse/fine flux-register accumulation
- [x] owner-local reflux, average-down, and temperature recovery
- [x] cross-owner adjacent PPM time-integrated face reconciliation
- [x] exact local and per-level cumulative subcycle accounting
- [x] serial parity for four-level branched and adjacent PPM trees
- [x] exact sparse rollback after a deep owner hydro failure
- [x] 1/2/4/8-rank Release and Debug gates

Chemistry and hydro now run directly on sparse AMR storage. Molecular
transport and combined full physics still use the replicated bridge.
Point-to-point communication and topology-changing distributed regrid remain
pending.

## Direct molecular transport on sparse MPI AMR patches (`0.57.0`)

- [x] owner-only SSPRK2 transport on rank-local sparse payloads
- [x] cumulative `r^2` recursive subcycling without a full tree replica
- [x] streamed parent interval states and effective diffusive fluxes
- [x] replicated compact diffusive flux-register accumulation
- [x] owner-local diffusive reflux, average-down, and temperature recovery
- [x] cross-owner adjacent shared diffusive face reconciliation
- [x] exact local and per-level parabolic subcycle accounting
- [x] serial parity for four-level branched and adjacent-child trees
- [x] exact sparse rollback after a deep owner transport failure
- [x] 1/2/4/8-rank Release and Debug gates

Chemistry, hydro, and molecular transport now all consume sparse AMR storage
directly. Sparse full-physics composition, point-to-point communication, and
topology-changing distributed regrid remain pending.

## Transactional full physics on sparse MPI AMR patches (`0.58.0`)

- [x] direct sparse `R-T-H-T-R` composition without a replicated-tree bridge
- [x] collective preflight for interval, hierarchy, and transport database
- [x] one outer sparse backup covering every operator stage
- [x] exact local and global chemistry, hydro, and transport call accounting
- [x] synchronized time, step, hydro, and transport subcycle counters
- [x] serial full-field parity and composite conservation on four levels
- [x] missing-transport rejection before mutation
- [x] exact outer rollback after chemistry and transport precede hydro failure
- [x] 1/2/4/8-rank Release and Debug gates

The complete qualified physics interval now consumes rank-local AMR payloads
directly. Component transactions retain their own backups, and their interval
data and compact registers still use collective communication. Reducing that
temporary memory and communication cost, point-to-point schedules, and
topology-changing distributed regrid remain pending.

## Topology-changing sparse MPI AMR regrid (`0.59.0`)

- [x] explicit-plan sparse regrid with collective topology agreement
- [x] unchanged-plan no-op with synchronized evaluation accounting
- [x] conservative average-down, prolongation, and exact fine-overlap transfer
- [x] deterministic owner-map rebuild for changed patch counts and geometry
- [x] one-copy global sparse storage after the topology transition
- [x] serial full-field parity and composite conservation
- [x] transactional sparse solution and distribution rollback on invalid plans
- [x] 1/2/4/8-rank Release and Debug gates

Topology-changing explicit plans can now be committed through the sparse MPI
API. The transition temporarily materializes a full correctness replica on
each rank before the rebuilt tree is scattered to its new owners. Direct
tag-driven sparse planning, point-to-point overlap transfer, and removal of
that temporary replica remain pending.

## Tag-driven sparse MPI AMR regrid (`0.60.0`)

- [x] solution-driven parent-local tag planning through four levels
- [x] disconnected-feature clustering with deterministic parent relationships
- [x] collective tagged-cell, topology-change, and overlap-count agreement
- [x] shared transactional topology commit and rebuilt owner distribution
- [x] unchanged tag-derived plan no-op with evaluation accounting
- [x] globally single-copy sparse storage after tag-driven rebuild
- [x] serial full-field parity and composite conservation
- [x] invalid-tag rollback for sparse solution and owner distribution
- [x] 1/2/4/8-rank Release and Debug gates

The sparse public API now covers both explicit and solution-derived topology
changes. Tag planning and conservative rebuild still operate on a temporary
all-rank correctness replica. Owner-local tagging, point-to-point overlap
transfer, and removal of that replica remain pending.

## Point-to-point sparse owner migration (`0.61.0`)

- [x] one packed old-owner to new-owner message per changed patch
- [x] exact state, temperature, narrow-ghost, and wide-ghost reconstruction
- [x] local copy without MPI traffic when ownership is unchanged
- [x] no payload allocation on ranks unrelated to a migrated patch
- [x] collective acceptance after every ordered direct transfer
- [x] global transfer count equals the changed-owner patch count
- [x] one-copy storage and exact replicated gather after migration
- [x] 1/2/4/8-rank Release and Debug gates

Same-hierarchy rebalancing now uses point-to-point patch payload transfer
instead of per-patch broadcast. Physics interval data, halo exchange, regrid
overlap transfer, and temporary topology replicas remain on the collective
correctness schedule.

## Point-to-point sparse adjacent halos (`0.62.0`)

- [x] one packed bidirectional exchange per cross-owner adjacent sibling face
- [x] one-layer narrow and four-layer PPM state/temperature payloads
- [x] local boundary copy without MPI traffic for same-owner siblings
- [x] no payload allocation on ranks unrelated to the adjacent face
- [x] exact global transfer count from the sparse chemistry refresh
- [x] unchanged cross-owner chemistry, hydro, and transport serial parity
- [x] unchanged composite conservation and transactional rollback gates
- [x] 1/2/4/8-rank Release and Debug gates

Adjacent same-level ghost traffic is now point to point on sparse storage.
Parent interval streaming, child-to-parent synchronization, flux
reconciliation, average-down, topology overlap transfer, and the temporary
regrid correctness replica remain collective work.

## Point-to-point sparse child-to-parent transfer (`0.63.0`)

- [x] one direct interior-state send per cross-owner child/parent pair
- [x] shared transfer path for chemistry average-down and physics synchronization
- [x] local interior copy without MPI traffic for same-owner pairs
- [x] no child payload allocation on unrelated ranks
- [x] exact global chemistry transfer count from the owner map
- [x] unchanged sparse chemistry, hydro, and transport serial parity
- [x] unchanged conservation, subcycle accounting, and rollback gates
- [x] 1/2/4/8-rank Release and Debug gates

Sparse child interiors now reach only the parent owner that consumes them.
Parent interval states and fluxes, parent-to-child ghost fill, shared-flux
reconciliation, topology overlap transfer, and the temporary regrid correctness
replica remain collective work.

## Point-to-point sparse parent-state fanout (`0.64.0`)

- [x] one direct parent-state send per distinct remote child owner
- [x] one received state reused for every local child of the same parent
- [x] same-owner parent state consumed without MPI traffic
- [x] no parent-state allocation on ranks unrelated to the parent or children
- [x] exact global recipient count derived independently from the owner map
- [x] unchanged sparse chemistry, hydro, and transport serial parity
- [x] unchanged conservation, subcycle accounting, and rollback gates
- [x] 1/2/4/8-rank Release and Debug gates

Final sparse parent-to-child ghost refresh is now point to point. Recursive
hydro/transport interval states, face fluxes, level counters, shared-flux
reconciliation, topology overlap transfer, and the temporary regrid correctness
replica remain collective work.

## Broadcast-free sparse recursive physics (`0.65.0`)

- [x] packed interval start/end state fanout to distinct child owners only
- [x] direct child boundary-flux return to the parent owner
- [x] parent-owner shared-flux construction and direct child correction
- [x] owner-local coarse/fine flux-register accumulation and reflux
- [x] one stage-boundary reduction for owner-local level-counter deltas
- [x] zero `MPI_Bcast` calls in the sparse physics module
- [x] exact communication counts for hyperbolic and parabolic subcycling
- [x] unchanged serial parity, conservation, call accounting, and rollback
- [x] 1/2/4/8-rank Release and Debug gates

Sparse chemistry, hydro, transport, full-physics composition, owner migration,
same-level halos, parent/child synchronization, and final ghost refresh now
operate without broadcast-based payload replication. Topology-changing regrid
still materializes and rebuilds through a temporary all-rank correctness
replica; direct distributed tagging and overlap transfer remain pending.

## Direct explicit-plan sparse regrid (`0.66.0`)

- [x] candidate hierarchy construction without materialized patch fields
- [x] owner-local parent-to-child prolongation through arbitrary depth
- [x] one direct state message per cross-owner new child
- [x] direct old-owner to new-owner fine-overlap segment transfer
- [x] state and temperature preservation for every retained overlap cell
- [x] exact independent prolongation and overlap message counts
- [x] unchanged-plan metadata-only no-op
- [x] exact serial full-field parity and composite conservation
- [x] transactional solution and owner-map rollback for invalid plans
- [x] 1/2/4/8-rank Release and Debug gates

Explicit-plan topology changes now remain sparse throughout hierarchy rebuild,
prolongation, overlap retention, average-down, and ghost refresh. No rank
materializes the complete old or new field tree. Tag-driven planning still uses
the temporary correctness replica and remains the next distributed regrid gap.

## Owner-local sparse tag planning (`0.67.0`)

- [x] gradient tagging only on each candidate parent owner
- [x] compact integer plan metadata agreement without field replication
- [x] owner-local candidate hierarchy prolongation through configured depth
- [x] exact tag-evaluation and candidate-transfer accounting
- [x] shared direct explicit-plan commit for the final tagged hierarchy
- [x] deterministic temperature recovery after conserved-field installation
- [x] exact tagged hierarchy, field, ghost, and bookkeeping serial parity
- [x] unchanged-plan no-op and invalid-tag transactional rollback
- [x] no `MPI_Bcast` or materialized-tree path in the sparse MPI module
- [x] 1/2/4/8-rank Release and Debug gates

Both explicit and solution-tagged sparse topology changes now keep field data
globally single-copy from planning through commit. Only compact hierarchy and
owner-map metadata are replicated; collective logical acceptance remains the
transaction boundary.

## General-EOS PeleC-style reactive Riemann solver (`0.68.0`)

- [x] left/right NASA7 frozen sound speed and acoustic impedance
- [x] PeleC star pressure and normal velocity estimate
- [x] upwind or stationary-averaged species-density origin state
- [x] EOS-checked star-density correction and inward/outward wave interpolation
- [x] final interface internal energy rebuilt from the NASA7 mixture EOS
- [x] exact species-flux closure to the total mass flux
- [x] equal-state, stationary/moving material-contact, and shock gates
- [x] x/y directional momentum rotation
- [x] invalid-state rejection and cleared optional interface outputs
- [x] 40/80/160-cell composition-wave second-order convergence
- [x] selectable `pelec` path in reactive 1D and 2D configurations
- [x] complete Release and Debug CI gates

The reactive acoustic path now follows PeleC `Source/Riemann.H` using the
available NASA7 ideal-gas-mixture EOS instead of the earlier constant-`gamma`
reduction. HLLC and Rusanov remain independent selectable baselines. Rotating
frames, auxiliary/linear advected fields, embedded boundaries, and arbitrary
PelePhysics EOS models remain outside this qualification.

## Prescribed species wall flux (`0.69.0`)

- [x] per-face impermeable/prescribed species-mode configuration
- [x] wall-to-gas input convention and lower/upper coordinate orientation
- [x] finite, active-species, and exact zero-net-mass validation
- [x] rejection on non-wall faces or without enabled species transport
- [x] species-enthalpy contribution to total-energy flux
- [x] shared positivity limiting of species and species enthalpy
- [x] x/y face unit gates, invalid-vector rejection, and application parsing
- [x] transient inventory, total-mass, cellwise closure, and temperature gates
- [x] complete Release and Debug CI gates

This is the transport boundary needed by later catalytic-wall models. It does
not yet calculate surface reaction rates, coverages, or wall chemistry.

## Subcycle-weighted MPI AMR ownership (`0.70.0`)

- [x] 64-bit per-patch and per-rank estimated work metadata
- [x] cell/storage exponent 0 compatibility mode
- [x] cumulative hyperbolic `r`-subcycle exponent 1
- [x] cumulative parabolic `r^2`-subcycle exponent 2
- [x] deterministic least-work owner selection with lowest-rank tie break
- [x] collective work-model agreement and invalid-exponent rejection
- [x] distribution validity checks for cell, patch, and work totals
- [x] explicit and owner-local tag-driven sparse regrid model preservation
- [x] reduced four-level maximum estimated work at 2 and 4 ranks
- [x] non-increasing maximum work at 1 and 8 ranks
- [x] unchanged sparse physics, migration, conservation, and rollback gates
- [x] 1/2/4/8-rank Release and Debug MPI gates

The model accounts for deterministic level subcycling, not measured wall-clock
time. Runtime feedback, heterogeneous CPU/GPU weights, and dynamic work
stealing remain outside this milestone.

## Distributed sparse MPI AMR timestep (`0.71.0`)

- [x] owner-local hyperbolic CFL evaluation on sparse patch payloads
- [x] optional owner-local molecular-transport stability evaluation
- [x] cumulative `r` hyperbolic conversion to the root-step limit
- [x] cumulative `r^2` parabolic conversion to the root-step limit
- [x] one communicator-wide minimum reduction without field gathering
- [x] neutral participation by ranks with no locally owned patches
- [x] collective invalid-state and missing-transport rejection
- [x] deterministic zero timestep on every rejected path
- [x] exact serial patch-tree timestep parity
- [x] 1/2/4/8-rank Release and Debug MPI gates

This supplies the distributed stability decision needed by a runnable sparse
MPI AMR driver. Adaptive stepping, stop-time clipping, regrid cadence, output,
and restart orchestration remain driver responsibilities.

## Runnable sparse MPI AMR application (`0.72.0`)

- [x] public namelist-driven `pelef_mpi_amr_reactive_1d` executable
- [x] configurable cell, hyperbolic, or parabolic ownership weighting
- [x] initial root scatter followed by owner-local arbitrary-depth tagging
- [x] distributed stop-time-clipped timestep and `R-T-H-T-R` loop
- [x] configured periodic owner-local topology rebuilds
- [x] final sparse-to-replicated diagnostic gather
- [x] physically ordered arbitrary-depth composite patch-tree CSV output
- [x] positive density, pressure, temperature, and species-closure checks
- [x] exact output parity at 1, 2, 4, and 8 ranks
- [x] GNU Fortran Release and bounds/FPE-checked Debug gates

The normal timestep and regrid path retains globally single-copy patch fields.
The initial root setup and final diagnostics/output materialize field data.

## Rank-independent sparse MPI AMR restart (`0.73.0`)

- [x] versioned, self-describing reactive patch-tree checkpoint format
- [x] mechanism species-name and conserved-width compatibility checks
- [x] base geometry, arbitrary-depth hierarchy, state, and temperature restore
- [x] physical time, step, regrid, overlap, and level-advance accounting restore
- [x] configurable checkpoint cadence and clean stop-after-write mode
- [x] owner-map-free restart with deterministic redistribution to active ranks
- [x] serial checkpoint round-trip regression
- [x] two-rank checkpoint resumed on four and eight ranks
- [x] uninterrupted/restarted composite field parity within `5e-13`
- [x] GNU Fortran Release and bounds/FPE-checked Debug gates

Scheduled checkpoints materialize the sparse tree and rank zero writes one
formatted file. Parallel plotfiles, atomic multi-file commits, asynchronous
I/O, and runtime-measured load balancing remain outside this milestone.

## Embedded-boundary geometry foundation (`0.74.0`)

- [x] nodal Cartesian level-set input with an explicit positive-fluid contract
- [x] bounded cell-volume fractions from clipped affine triangles
- [x] shared x/y face open-area fractions from endpoint interpolation
- [x] regular, cut, and covered cell classification
- [x] self-validating extents, spacing, array bounds, finite values, and types
- [x] exact all-fluid and all-solid geometry gates
- [x] machine-precision vertical and diagonal planar-interface area gates
- [x] refined circular-interface integrated-area convergence gate
- [x] GNU Fortran Release and bounds/FPE-checked Debug gates

Cut-face centroids/normals, cut-cell flux divergence, state redistribution,
small-cell stabilization, embedded-wall conditions, AMR coupling, and MPI
distribution remain separate milestones.

## Embedded-boundary interface metrics (`0.75.0`)

- [x] physical embedded-boundary length in every cut cell
- [x] length-weighted physical interface centroid
- [x] unit `grad(phi)` normal directed from solid toward fluid
- [x] duplicate suppression for interfaces coincident with the cell diagonal
- [x] exact vertical and diagonal planar length/centroid/normal gates
- [x] circular perimeter error reduction from `20x20` to `40x40`
- [x] circular inward-fluid normal orientation and refinement gate
- [x] Release and bounds/FPE-checked Debug qualification

These metrics enable the next cut-wall flux milestone. Flux divergence,
small-cell stabilization, wall thermodynamics, AMR, and MPI remain unclaimed.

## Reactive embedded-boundary slip-wall flux (`0.76.0`)

- [x] general-EOS pressure recovery from the reactive conserved state
- [x] arbitrary-orientation stationary slip-wall momentum flux
- [x] explicit solid-to-fluid versus fluid-outward normal sign contract
- [x] cut-interface-length integration and fluid-volume normalization
- [x] exact zero mass, tangential-z momentum, energy, and species wall fluxes
- [x] rotational covariance and velocity-independence gates
- [x] vertical and diagonal integrated pressure-force balance
- [x] rejected-input zero-output transaction contract
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

This is the embedded-wall contribution to a cut-cell finite-volume operator.
Cartesian open-face divergence, a complete state update, small-cell
stabilization, viscous/thermal/catalytic walls, AMR, and MPI remain separate.

## Conservative reactive EB flux divergence (`0.77.0`)

- [x] solid-to-fluid integrated interface-normal vector in every cut cell
- [x] separate unit-normal and integrated-normal geometry contracts
- [x] shared Cartesian fluxes weighted by x/y open-face fractions
- [x] wall and Cartesian contributions divided by the same fluid volume
- [x] exactly inert covered cells
- [x] regular, vertical, diagonal, and circular uniform-pressure balance gates
- [x] finite face-flux and exact array-extent validation
- [x] rejected-input zero-output transaction contract
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

This supplies a conservative semidiscrete cut-cell operator. Flux generation
near covered neighbors, a complete state update, small-cell stabilization,
viscous/thermal/catalytic walls, AMR, and MPI remain separate milestones.

## Conservative EB flux redistribution and state update (`0.78.0`)

- [x] face-connected active-cell redistribution neighborhoods
- [x] volume-fraction-weighted nonconservative neighborhood update
- [x] stable cut-cell blend with extensive excess-update redistribution
- [x] exact volume-weighted conservation for every conserved component
- [x] exact preservation of uniform active-cell right-hand sides
- [x] exactly inert covered cells and compact one-face support
- [x] EOS-validated reactive forward update with full rollback
- [x] positive update at volume fraction `0.05` where the raw update is negative
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

This qualifies the documented first-order AMReX FluxRedist construction on a
single Cartesian level. PeleC's default weighted StateRedist, boundary/periodic
redistribution neighborhoods, higher-order reconstruction, EB flux generation,
AMR coupling, and MPI distribution remain separate milestones.

## Weighted EB state redistribution (`0.79.0`)

- [x] AMReX-style aperture-normal merge-neighborhood construction
- [x] default target volume fraction `0.5`
- [x] one, two, and diagonal third-neighbor selection in two dimensions
- [x] explicit overlapping-neighborhood count (`nrs`)
- [x] partitioned self and neighbor weights (`alpha`)
- [x] zeroth-order weighted neighborhood states
- [x] exact uniform-state preservation and componentwise conservation gates
- [x] shared-cell analytical gate for two overlapping small-cell neighborhoods
- [x] EOS-validated provisional reactive-state update with full rollback
- [x] positive update at volume fraction `0.05` where the raw state is negative
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

This matches the current PeleC default combination of `StateRedist` and
`eb_srd_max_order = 0` on one nonperiodic Cartesian level. Higher-order
reconstruction and limiting, periodic/physical ghost-cell neighborhoods,
multilevel redistribution, EB flux generation, AMR coupling, and MPI
distribution remain separate milestones.

## Complete first-order reactive EB hydro update (`0.80.0`)

- [x] reactive Riemann fluxes only on positive-aperture Cartesian faces
- [x] exactly zero closed-face fluxes
- [x] zero-gradient lower and upper domain-face states
- [x] rejection of a positive-aperture face touching a covered cell
- [x] open-area and integrated slip-wall conservative divergence
- [x] PeleC-default weighted StateRedist on the provisional state
- [x] general-EOS recovery and whole-step rollback
- [x] direct interior-face HLLC parity gate
- [x] regular, vertical, diagonal, and circular stationary-pressure gates
- [x] unknown-solver and nonfinite-state transaction gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

This connects the previously independent EB kernels into one first-order
reactive Euler update with nonperiodic zero-gradient outer boundaries. It does
not yet claim high-order cut-cell reconstruction and transverse prediction,
face-centroid flux interpolation, periodic boundaries, moving walls, transport
or chemistry splitting, AMR coupling, or MPI distribution.

## Runnable reactive EB 2D application (`0.81.0`)

- [x] combined `reactive_2d` and `embedded_boundary` namelist input
- [x] configurable plane and inside/outside-circle level sets
- [x] explicit rejection of geometry without cut cells
- [x] active-cell general-EOS CFL timestep with exact final-time clipping
- [x] complete repeated PCM/HLLC/StateRedist hydro advancement
- [x] componentwise volume-fraction-weighted conserved integrals
- [x] active-cell density, pressure, temperature, speed, and closure extrema
- [x] CSV cell centers, volume fraction/type, wall metrics, and flow fields
- [x] committed circular-obstacle input and output-structure regression
- [x] direct-API rejection of unsupported chemistry and transport modes
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The standalone `pelef_reactive_eb_2d` path advances stationary embedded
boundaries through a real input-to-output workflow. It is deliberately limited
to first-order hydro with outflow outer boundaries. Chemistry, molecular
transport, higher-order reconstruction, AMR, and MPI remain unclaimed.

## Active-cell EB chemistry splitting (`0.82.0`)

- [x] optional active mask on the shared reactive 2D chemistry operator
- [x] candidate-array chemistry commit with no partial-cell mutation
- [x] covered cells excluded from primitive recovery and reactor integration
- [x] reaction--EB hydro--reaction Strang composition
- [x] outer rollback when hydro fails after the first reaction half-step
- [x] elementary and full-mechanism dispatch in the public EB application
- [x] active-cell field parity with the regular 2D chemistry application
- [x] bitwise covered-cell parity with a chemistry-disabled EB run
- [x] volume-weighted mass and total-energy conservation gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The EB application is now genuinely reactive for stationary adiabatic slip
walls. Molecular transport, higher-order EB hydro, AMR coupling, and MPI
ownership remain outside this serial single-level path.

## Active-stencil EB PLM and face-centroid fluxes (`0.83.0`)

- [x] normalized tangential centroid offset for every Cartesian EB face
- [x] PCM-compatible face-center Godunov flux staging
- [x] frozen-composition characteristic PLM normal prediction
- [x] two-sided active-cell slope contract with local zero-slope fallback
- [x] density, pressure, and composition slope scaling and sanitization
- [x] linear face-center to open-face-centroid flux interpolation
- [x] exact affine moving-contact mass-flux gates in x and y
- [x] regular, vertical, diagonal, and circular uniform-state PCM/PLM parity
- [x] face-centroid sign, magnitude, interpolation, and validity gates
- [x] unknown-limiter zero-output and complete-step rollback gates
- [x] input-driven characteristic-PLM circular-obstacle regression
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

This adds the PeleC-ordered face-center Godunov then tangential face-centroid
interpolation path while keeping PCM as the default API behavior. It does not
claim PeleC's unsplit EB transverse predictor, PPM, higher-order StateRedist,
periodic EB neighborhoods, molecular transport, AMR, or MPI ownership.

## Second-order weighted EB StateRedist (`0.84.0`)

- [x] normalized fluid-volume centroids from the clipped positive polygons
- [x] partition-consistent weighted neighborhood centroids (`cent_hat`)
- [x] selectable `state_redist_max_order = 0` or `2`
- [x] active connected 3-by-3 least-squares `Qhat` slopes
- [x] active grown 5-by-5 fallback for a rank-deficient normal matrix
- [x] AMReX-style pairwise centroid limiter
- [x] merge-recipient maximum-principle slope scale
- [x] zero-moment linear scatter through overlapping neighborhoods
- [x] exact affine-state reproduction and volume-weighted conservation gates
- [x] discontinuous-state bounds, invalid-order, and rollback gates
- [x] input-driven second-order circular and chemistry EB regressions
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The API keeps zeroth-order StateRedist as its compatibility default, while the
committed EB applications select second order. This milestone does not claim
fourth-order slopes, physical/periodic ghost neighborhoods, multilevel
redistribution, EB coarse/fine flux synchronization, molecular transport, or
MPI ownership.

## Static two-level EB average-down (`0.85.0`)

- [x] aligned rectangular fine-patch metadata and geometry validation
- [x] parent/child EB volume-measure consistency gate
- [x] AMReX-style fine-volume-fraction-weighted state restriction
- [x] first-child fallback for a fully covered fine block
- [x] outside-patch coarse-state preservation
- [x] composite EB integral with single-count coarse/fine coverage
- [x] exact constant and affine fluid-centroid restriction gates
- [x] composite-to-restricted conservation gate
- [x] reactive active-parent EOS and temperature recovery
- [x] covered-parent preservation and whole-array rollback
- [x] nonfinite, nonphysical, and misaligned-hierarchy rejection gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

This is a serial synchronization foundation for one static two-level rectangle.
It does not yet advance an EB hierarchy or provide prolongation, ghost fill,
subcycling, coarse/fine flux registers and reflux, multiple fine patches,
dynamic regridding, multilevel redistribution, or MPI ownership.

## EB coarse/fine flux register and cut-cell re-reflux (`0.86.0`)

- [x] coarse/fine boundary face-measure compatibility gate
- [x] independent time-weighted coarse and fine flux accumulation
- [x] four-sided open-area and physical-subface flux matching
- [x] exact one-coarse/two-fine subcycle cancellation gate
- [x] regular exterior-cell reflux
- [x] AMReX-style cut-cell `kappa` re-reflux stabilization
- [x] connected 3-by-3 fluid-volume recipient distribution
- [x] fine-child transfer for recipients below the fine rectangle
- [x] raw-to-redistributed composite conservation gate
- [x] register reset only after successful state commit
- [x] two-level reactive state and temperature EOS transaction
- [x] covered-state preservation and nonphysical rollback gate
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

This completes the conservative synchronization kernels for one static serial
two-level EB rectangle. It does not yet build or advance the hierarchy, fill
coarse/fine ghosts, prolong initial data, drive fine substeps, regrid, support
multiple patches or deeper levels, or distribute EB ownership with MPI.

## Static two-level reactive EB hydrodynamic advance (`0.87.0`)

- [x] piecewise-constant coarse-to-fine state initialization
- [x] active fine-child EOS temperature recovery
- [x] one coarse step paired with `r` fine hyperbolic substeps
- [x] conserved coarse-time interpolation at fine-patch exterior faces
- [x] EOS recovery for every open exterior face state
- [x] PCM or characteristic-PLM fine-level reconstruction
- [x] configurable weighted-StateRedist target volume fraction
- [x] exact advancing centroid fluxes supplied to the EB register
- [x] reactive cut-cell re-reflux followed by EB average-down
- [x] hierarchy-wide state and temperature transaction
- [x] constant-state prolong/restrict and subcycled preservation gates
- [x] nonmatching interface mass, energy, and species conservation gate
- [x] composite conservation and covered-state preservation gates
- [x] invalid interpolation-time and solver rollback gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

This is the first complete hydrodynamic interval for one static, aligned,
strictly internal fine rectangle. Coarse exterior data is piecewise constant in
space and interpolated in time. Chemistry, molecular transport, automatic
timestep selection, dynamic regridding, multiple patches, deeper levels, and
MPI ownership are not yet composed with this EB hierarchy.

## Runnable static two-level reactive EB AMR application (`0.88.0`)

- [x] `eb_amr` namelist for internal patch bounds and refinement ratio
- [x] shared configured plane/circle level-set builder on arbitrary regions
- [x] independently reconstructed coarse and fine EB geometry
- [x] qualified patch construction and PCM fine initialization
- [x] coarse stability limit from both coarse and fine active states
- [x] final-time clipping and maximum-step rejection
- [x] repeated transactional two-level hydrodynamic advance
- [x] initial/final composite conserved diagnostics
- [x] separate synchronized coarse and fine geometry/state CSV output
- [x] installed `pelef_reactive_eb_amr_2d` executable
- [x] input-driven diagonal-EB regression and CSV structural checker
- [x] direct unsupported-chemistry rejection gate
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The static hierarchy is now runnable rather than only callable as a kernel.
The app remains hydrodynamics-only and owns one strictly internal aligned fine
rectangle. Chemistry, molecular transport, automatic tagging/regridding,
multiple patches, deeper levels, checkpoint/restart, and MPI ownership remain
future EB AMR work.

## Solution-driven two-level reactive EB AMR regridding (`0.89.0`)

- [x] active-cell four-neighbor temperature-gradient tagging
- [x] combined relative threshold and absolute noise floor
- [x] covered-cell and physical-boundary exclusion
- [x] buffered, minimum-size, strictly internal rectangular patch planning
- [x] old fine-patch reactive EB average-down before topology replacement
- [x] new fine-patch PCM initialization from the synchronized root
- [x] exact same-resolution fine-state and temperature overlap retention
- [x] active new-fine EOS validation and whole-hierarchy transaction
- [x] initial and periodic accepted-step regrid cadence
- [x] unchanged-patch and empty-tag retention without false regrid counts
- [x] public namelist controls and committed regrid count diagnostics
- [x] moving-hotspot application regression and geometry-aware CSV checker
- [x] composite conservation, new-cell PCM, retired-cell restriction, and
  nonfinite rollback gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The serial application now moves or resizes one ratio-aligned fine rectangle
around root temperature structure. It remains hydrodynamics-only and retains
the existing patch when tags disappear. Multiple patches, deeper levels,
fine-patch removal, EB AMR chemistry/transport, checkpoint/restart, and MPI
ownership remain future work.

## Optional reactive EB AMR fine-patch lifecycle (`0.90.0`)

- [x] optional initial regrid evaluation and untagged-patch removal policy
- [x] transactional fine-to-root collapse through reactive EB average-down
- [x] fine state, temperature, geometry, and patch-metadata release
- [x] root-only active-cell CFL selection and reactive EB hydro advance
- [x] PCM fine-patch re-creation when root temperature tags return
- [x] lifecycle-aware composite/root conserved diagnostics
- [x] inactive fine-output suppression and application diagnostics
- [x] time-loop degrid regression with released-storage verification
- [x] direct create-collapse conservation and EOS gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The public serial driver may now transition between one root level and one
strictly internal ratio-aligned fine rectangle. Multiple simultaneous patches,
deeper levels, EB AMR chemistry/transport, checkpoint/restart, and MPI ownership
remain future work.

## Reactive EB AMR chemistry composition (`0.91.0`)

- [x] elementary and full H2/O2 mechanism loading in the EB AMR application
- [x] active-cell reaction half-steps on both coarse and fine levels
- [x] `reaction-hydro-reaction` composition around EB subcycling and reflux
- [x] post-reaction reactive fine-to-coarse average-down
- [x] covered-cell chemistry exclusion on independently built EB geometries
- [x] hierarchy-wide state and temperature rollback after a later hydro failure
- [x] root-only chemistry through the qualified single-level EB Strang path
- [x] reacting fine-patch-collapse lifecycle with mass and energy gates
- [x] regular-grid reference parity on active coarse and fine EB AMR cells
- [x] explicit unsupported-molecular-transport rejection before any step
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The serial application now advances chemistry on either its active two-level
hierarchy or its root-only lifecycle state. The qualified topology remains one
strictly internal ratio-aligned fine rectangle. EB AMR molecular transport,
multiple simultaneous patches, deeper levels, checkpoint/restart, and MPI
ownership remain future work.

## Reactive EB AMR checkpoint/restart (`0.92.0`)

- [x] versioned formatted checkpoint schema and end marker
- [x] species-order and reactive-state-size compatibility checks
- [x] mesh, EB geometry, chemistry, hydro, redistribution, and regrid signature
- [x] coarse state plus optional active fine state and actual patch bounds
- [x] root-only lifecycle encoding without fine storage
- [x] time, step, regrid, minimum-timestep, and base-density metadata
- [x] geometry reconstruction and active-cell EOS temperature recovery
- [x] private-candidate read with malformed-file transactional rejection
- [x] periodic/final writes and optional stop-after-checkpoint control
- [x] input-driven restart with mutable final time, step budget, and outputs
- [x] uninterrupted versus reacting fine-to-root restart field parity
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

Longer serial EB AMR runs can now stop after a committed checkpoint and resume
without recreating an obsolete fine patch or losing lifecycle counters. The
schema is intentionally serial and formatted; distributed checkpoint I/O,
multiple simultaneous patches, deeper levels, and EB AMR molecular transport
remain future work.

## Reactive EB AMR multipatch kernel (`0.93.0`)

- [x] deterministic disconnected-tag clustering with configurable gap joining
- [x] per-cluster buffer/minimum-size expansion and strict boundary rejection
- [x] candidate coalescing for the 3-by-3 redistribution separation contract
- [x] ordered multipatch geometry, state, temperature, and validity ownership
- [x] conservative PCM creation and exact old/new fine-overlap retention
- [x] conservative patch-set movement, repartition, removal, and average-down
- [x] composite integration without double-counting multiply refined parents
- [x] one root advance and refinement-ratio subcycling for every child
- [x] independent EB flux registers, sequential reflux, and final average-down
- [x] hierarchy-wide hydrodynamic rollback with failure-stage diagnostics
- [x] active-cell Strang chemistry on the root and every child
- [x] post-reaction synchronization and whole-hierarchy rollback
- [x] mass, energy, species, overlap, synchronization, and failure gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The serial kernel now represents and advances multiple separated fine
rectangles over one EB root. The public EB AMR application, its runtime
configuration, CSV layout, and checkpoint schema still select the qualified
single-patch lifecycle. Connecting this patch set to that lifecycle, deeper
levels, EB AMR molecular transport, and MPI ownership remain separate work.

## Reactive EB AMR multipatch application (`0.94.0`)

- [x] `multipatch_enabled` runtime dispatch and maximum tag-gap input
- [x] configured seed rectangle converted to transactional patch-set ownership
- [x] initial and periodic disconnected-temperature-tag regridding
- [x] empty-set removal policy and unchanged-collection no-op behavior
- [x] root plus all-child active-cell CFL timestep selection
- [x] public set-wide Strang chemistry/hydro time loop and regrid counters
- [x] deterministic per-child CSV names, geometry, state, and diagnostics
- [x] explicit multipatch checkpoint/restart rejection before initialization
- [x] two-hotspot plane-EB input case producing two separated fine rectangles
- [x] public mass, energy, positivity, species-closure, and output gates
- [x] focused application gate before the complete CI regression
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The serial EB AMR executable can now select either the existing single-patch
checkpoint-capable lifecycle or a dynamic two-level patch set. Multipatch
checkpoint/restart, EB molecular transport, deeper levels, physical-boundary
fine patches, and distributed ownership remain separate work.

## Reactive EB AMR multipatch checkpoint/restart (`0.95.0`)

- [x] dedicated versioned patch-set checkpoint magic and schema
- [x] species, mesh, EB, physics, regrid, and collection compatibility signature
- [x] exact ordered child count, bounds, dimensions, state, and temperature
- [x] root state, time, step/regrid counters, minimum timestep, and base density
- [x] private geometry rebuild and EOS temperature recovery for every level
- [x] complete-set separation and terminal-marker validation before publication
- [x] scheduled post-regrid writes, stop-after-write, and cadence-preserving restart
- [x] direct two-child round trip and truncated-file rollback gates
- [x] public reacting reference, checkpoint-stop, restart, and CSV parity gate
- [x] unchanged single-patch checkpoint format and restart behavior
- [x] focused checkpoint gate before the complete CI regression
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The public serial EB AMR executable can now stop and resume either the default
single-patch lifecycle or an explicitly enabled patch set without losing
accepted topology or cadence. EB molecular transport, deeper levels,
physical-boundary fine patches, and distributed ownership remain separate work.

## Reactive EB AMR configured physical-boundary patch (`0.96.0`)

- [x] domain-inclusive configured single-patch input and geometry construction
- [x] outflow-side exterior state copied from the current fine boundary cell
- [x] coarse-time interpolation retained on every coarse/fine patch side
- [x] physical-side flux-register accumulation and reflux omitted
- [x] fine exterior payload dimension, finiteness, and temperature validation
- [x] public x-lower-boundary plane-EB application case
- [x] exact coarse/fine coordinates, final time, and EB-class output gates
- [x] finite positive thermodynamics, stationary uniform state, and closure gate
- [x] focused application gate before the complete CI regression
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The public serial EB AMR executable can now advance a configured static fine
rectangle that meets an outflow physical boundary. Dynamic tag planning still
requires strictly internal patches. EB molecular transport, deeper EB levels,
physical-side dynamic clustering, non-outflow refined boundaries, and
distributed ownership remain separate work.

## Reactive EB AMR dynamic physical-boundary planning (`0.97.0`)

- [x] one-sided temperature-gradient tags on active root boundary cells
- [x] domain-inclusive single-patch plan validation, buffering, and growth
- [x] domain-inclusive multipatch flood fill, clustering, and coalescing
- [x] full-domain minimum patch sizes accepted by input and runtime validation
- [x] physical-side child support in patch-set hydro and Strang transactions
- [x] conservative topology movement away from a physical side
- [x] unit gates for boundary tags, plans, multipatch hydro, and overlap transfer
- [x] public boundary-plus-interior double-hotspot dynamic application
- [x] aligned child output, separation, EB classes, positivity, and closure gates
- [x] focused application gate before the complete CI regression
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The serial EB AMR application can now dynamically refine an outflow boundary
with either one child or a separated patch set. EB molecular transport, deeper
EB levels, non-outflow refined boundaries, distributed ownership, and explicit
physical-side checkpoint parity remain separate work.

## Static three-level EB synchronization (`0.98.0`)

- [x] strictly nested root, middle, and finest EB patch validation
- [x] three-level composite integral with exact finest-level ownership
- [x] deepest-to-middle then middle-to-root volume-weighted average-down
- [x] generic conserved-field transaction with no partial publication
- [x] reactive state restriction and EOS temperature recovery at both parents
- [x] unchanged uncovered root and middle ownership
- [x] nonfinite-state and invalid-temperature rollback gates
- [x] composite conservation across regular, cut, and covered EB cells
- [x] focused gate before the complete CI regression
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The EB AMR transfer layer now synchronizes one static three-level hierarchy.
It does not yet own a three-level lifecycle or perform recursive hydro
subcycling, flux-register reflux, chemistry, regridding, checkpointing, output,
transport, or MPI distribution.

## Static three-level reactive EB hydrodynamics (`0.99.0`)

- [x] one root, `r1` middle, and `r1*r2` finest hydro updates
- [x] parent-time interpolation at both nested interfaces
- [x] independent root/middle and middle/finest EB flux registers
- [x] inner reflux and synchronization after every middle interval
- [x] outer reflux followed by deepest-to-root EOS synchronization
- [x] composite mass, total-energy, and species conservation gate
- [x] finite positive temperatures and final synchronization gate
- [x] whole-hierarchy solver-failure rollback
- [x] two-cell finest separation and regular-interface validation
- [x] explicit EB-cut finest-interface rejection
- [x] focused gate before the complete CI regression
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The static serial kernel now advances a three-level reactive EB hierarchy when
the finest coarse/fine interface lies entirely in regular fluid. EB-cut nested
interfaces require a dedicated multilevel geometric reflux correction and are
rejected. Chemistry composition, a public lifecycle, regridding,
checkpoint/restart, output, transport, arbitrary depth, and MPI ownership
remain separate work.

## EB-cut nested-interface conservation closure (`0.100.0`)

- [x] pre-update middle/finest composite reference
- [x] signed time-integrated middle exterior flux target
- [x] post-reflux and post-average-down residual measurement
- [x] density, total-energy, and every-species correction
- [x] explicit exclusion of wall-force momentum components
- [x] density/species residual-closure validation
- [x] volume-weighted uncovered active-middle recipients
- [x] EOS temperature recovery for every corrected recipient
- [x] cut-interface three-level conservation and positivity gate
- [x] whole-hierarchy rollback on any closure or EOS failure
- [x] focused gate before the complete CI regression
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The static three-level hydro kernel now accepts a finest coarse/fine interface
crossed by the embedded boundary while retaining composite mass, total energy,
and species. The correction is conservative and transactional but globally
distributed over eligible middle cells; it is not yet PeleC-style local
multilevel EB redistribution. Chemistry, lifecycle ownership, regridding,
checkpointing, output, transport, arbitrary depth, and MPI remain separate.

## Static three-level reactive EB Strang chemistry (`0.101.0`)

- [x] active-cell reaction half-step on root, middle, and finest levels
- [x] recursive three-level EB hydro between reaction half-steps
- [x] compatibility with an EB-cut finest coarse/fine interface
- [x] post-chemistry finest-to-middle-to-root reactive average-down
- [x] composite mass and total-energy conservation
- [x] composite species activity with density/species closure
- [x] finite positive temperature on all three levels
- [x] whole-hierarchy rollback after chemistry when hydro rejects
- [x] focused gate before the complete CI regression
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The static serial driver now advances chemistry and recursively subcycled
hydrodynamics as one rollback-safe three-level EB transaction. It still does
not own a public three-level lifecycle, CFL loop, regridding, checkpointing,
output, molecular transport, arbitrary depth, or distributed MPI state.

## Public static three-level reactive EB AMR (`0.102.0`)

- [x] namelist-selected root/middle/finest static hierarchy
- [x] validated finest bounds with a two-middle-cell margin
- [x] root-to-middle and middle-to-finest PCM initialization
- [x] three-level active-cell CFL reduction with subcycle scaling
- [x] final-time clipping and accepted-step accounting
- [x] public recursive chemistry/hydro time loop
- [x] separate root, middle, and finest geometry-aware CSV output
- [x] cut-cell coverage and uniform-reactor reference parity on all levels
- [x] composite mass and total-energy conservation
- [x] focused application gate before the complete CI regression
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The public serial EB AMR executable can now run the qualified static
three-level hierarchy from a namelist. Dynamic three-level topology,
checkpoint/restart, molecular transport, arbitrary depth, and MPI ownership
remain separate work.

## Static three-level EB AMR checkpoint/restart (`0.103.0`)

- [x] dedicated versioned three-level checkpoint magic and schema
- [x] complete root, middle, and finest state and temperature payloads
- [x] nested topology, mechanism order, physics, and EB compatibility checks
- [x] private transactional read and all-level EOS temperature recovery
- [x] periodic/final writes and optional stop-after-write control
- [x] accepted time, minimum timestep, step count, and base-density recovery
- [x] uninterrupted versus restarted root/middle/finest CSV parity
- [x] unchanged single-patch and patch-set checkpoint schemas
- [x] focused application gate before the complete CI regression
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The public static three-level lifecycle can now stop after a committed
checkpoint and resume all three levels without reinitialization. Dynamic
three-level topology, molecular transport, arbitrary depth, distributed I/O,
and MPI ownership remain separate work.

## Tag-driven dynamic three-level finest patch (`0.104.0`)

- [x] active-middle-cell temperature-gradient tagging
- [x] two-cell-safe interior planning region
- [x] conservative old-finest average-down before topology replacement
- [x] PCM initialization of newly refined cells
- [x] exact state retention on overlapping old/new finest cells
- [x] transactional middle/finest publication and EOS validation
- [x] initialization-time and accepted-step regrid cadence
- [x] public committed-regrid accounting and failure-stage diagnostics
- [x] EB-aware 22 by 28 finest topology from an initial 8 by 8 patch
- [x] focused application gate before the complete CI regression
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The public three-level lifecycle can now move and resize its finest rectangle
inside a fixed middle level. Finest removal, dynamic middle/root topology,
sibling finest patches, checkpoint/restart after topology changes, molecular
transport, arbitrary depth, and MPI ownership remain separate work.

## Dynamic three-level topology checkpoint/restart (`0.105.0`)

- [x] distinct dynamic three-level checkpoint magic and schema
- [x] committed finest bounds and refinement-ratio persistence
- [x] accepted time, minimum timestep, step and regrid-count recovery
- [x] regrid cadence, tag threshold, buffer, and size compatibility checks
- [x] private reconstruction of the stored finest geometry
- [x] all-level EOS temperature recovery and transactional publication
- [x] unchanged static three-level, single-patch, and patch-set schemas
- [x] uninterrupted versus restarted root/middle/finest field parity
- [x] focused application gate before the complete CI regression
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

Dynamic three-level restart now resumes the committed finest topology rather
than the configured seed and retains regrid accounting across the split.
Finest removal or siblings, dynamic middle/root topology, molecular transport,
arbitrary depth, parallel checkpoint I/O, and MPI ownership remain separate
work.

## Single-level reactive EB molecular transport (`0.106.0`)

- [x] mixture viscosity, thermal conduction, and species diffusion flux reuse
- [x] barodiffusion and correction-velocity species mass closure
- [x] Cartesian-face to EB face-centroid interpolation
- [x] open-area and fluid-volume conservative divergence
- [x] adiabatic slip and species-impermeable embedded-wall closure
- [x] EB fluid-inventory species positivity limiter
- [x] SSPRK2 StateRedist and all-active-cell EOS recovery
- [x] transport stability limit in the public timestep selection
- [x] transactional chemistry/transport/hydro symmetric composition
- [x] public hotspot conduction and nonuniform conservation gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The public single-level EB application now runs the qualified molecular
transport subset. Thermal, viscous, catalytic, or prescribed-flux embedded
walls, EB AMR diffusive reflux, dynamic transport topology, MPI ownership, and
parallel transport remain separate work.

## Two-level reactive EB AMR molecular transport (`0.107.0`)

- [x] coarse and fine EB transport SSPRK2 transactions
- [x] ratio-subcycled fine transport intervals
- [x] time-interpolated coarse exterior states on fine patch boundaries
- [x] reusable EB face-centroid diffusive flux outputs
- [x] time-integrated coarse/fine diffusive flux register
- [x] conservative reflux and reactive average-down after each Euler stage
- [x] transport stability limits from both levels
- [x] transactional `R-T-H-T-R` driver composition
- [x] missing transport-database rejection and exact rollback
- [x] public inert/conducting hotspot comparison gate
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The public single-patch two-level EB AMR lifecycle now runs the qualified
mixture molecular-transport subset. Three-level and multipatch transport,
coarse-to-fine spatial slopes, non-outflow refined boundaries, thermal or
catalytic embedded walls, transport checkpoint/restart, MPI ownership, and
parallel transport remain separate work.

## Three-level reactive EB AMR molecular transport (`0.108.0`)

- [x] nested root/middle/finest EB transport SSPRK2 transactions
- [x] cumulative ratio subcycling on middle and finest levels
- [x] time-interpolated parent exterior states at both interfaces
- [x] independent time-integrated diffusive flux registers per interface
- [x] innermost reflux and EB-cut conservation closure per middle substep
- [x] outer reflux and deepest-first reactive average-down
- [x] root-equivalent parabolic stability limits from all three levels
- [x] transactional `R-T-H-T-R` three-level driver composition
- [x] missing transport-database rejection and exact rollback
- [x] public inert/conducting three-level hotspot comparison gate
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The public strictly nested three-level EB AMR lifecycle now runs the qualified
mixture molecular-transport subset for static or movable finest rectangles.
Multipatch transport, coarse-to-fine spatial slopes, non-outflow refined
boundaries, thermal or catalytic embedded walls, transport checkpoint/restart,
MPI ownership, and parallel transport remain separate work.

## Multipatch reactive EB AMR molecular transport (`0.109.0`)

- [x] one coarse EB transport update shared by all sibling patches
- [x] independent time-integrated diffusive flux register per child
- [x] ratio-subcycled fine transport with coarse-time exterior interpolation
- [x] sequential disjoint-interface reflux and set-wide reactive average-down
- [x] one global EB-cut composite conservation closure per Euler stage
- [x] density/species-consistent correction and EOS temperature recovery
- [x] SSPRK2 transaction with exact whole-patch-set rollback
- [x] all-child root-equivalent parabolic stability limits
- [x] transactional `R-T-H-T-R` multipatch driver composition
- [x] public inert/conducting double-hotspot comparison gate
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The public dynamic two-level sibling-patch EB AMR lifecycle now runs the
qualified mixture molecular-transport subset. Same-level diffusive exchange
for touching siblings, coarse-to-fine spatial slopes, non-outflow refined
boundaries, thermal or catalytic embedded walls, transport checkpoint/restart,
MPI ownership, and parallel transport remain separate work.

## MPI ownership foundation for reactive EB AMR 2D (`0.110.0`)

- [x] replicated root and sibling-patch topology consensus
- [x] contiguous nonoverlapping root y-tiles for every tested rank
- [x] unique deterministic owner for every root tile and fine child
- [x] raw, hyperbolic, or parabolic subcycle-weighted work models
- [x] 64-bit entity and per-rank work accounting
- [x] greedy deterministic load assignment with stable tie breaking
- [x] owner-authoritative root state and temperature synchronization
- [x] owner-authoritative child state and temperature synchronization
- [x] collective inconsistent-work-model rejection
- [x] transactional rollback for an invalid owner map
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The EB AMR MPI layer now owns distribution metadata and a correctness-first
replicated synchronization bridge. Direct rank-local EB physics, sparse field
storage, point-to-point coarse/fine traffic, distributed flux registers,
topology migration, parallel checkpoint I/O, and a public distributed time
loop remain separate work.

## Owner-only MPI reactive EB AMR chemistry (`0.111.0`)

- [x] replicated chemistry-control and mechanism-width consensus
- [x] active-cell root chemistry split across owned y-tiles
- [x] one active-cell chemistry update on each exclusive child owner
- [x] collective acceptance after every owner reactor transaction
- [x] owner-authoritative state and recovered-temperature broadcasts
- [x] replicated fine-to-root reactive average-down after chemistry
- [x] exact local and global committed entity-call accounting
- [x] serial patch-set chemistry and rank-count field parity
- [x] exact all-rank rollback after a late root-owner rejection
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

Reactive source integration now executes only on the owner of each root tile
or fine sibling, while the correctness bridge retains replicated field
storage. Owner-only EB hydrodynamics and molecular transport, sparse fields,
point-to-point coarse/fine traffic, distributed flux registers, topology
migration, parallel checkpoint I/O, and a public distributed time loop remain
separate work.

## Owner-only MPI reactive EB AMR hydrodynamics (`0.112.0`)

- [x] collective solver, reconstruction, limiter, timestep, and StateRedist consensus
- [x] one complete root-level EB finite-volume update on an exclusive physics owner
- [x] one ratio-subcycled fine update on each exclusive child owner
- [x] owner-local coarse/fine time-integrated flux registers
- [x] sequential cross-owner child reflux with authoritative broadcasts
- [x] root-owner fine-to-root reactive average-down
- [x] exact owner and global committed level-advance accounting
- [x] serial multipatch EB hydro field parity
- [x] exact all-rank rollback after a late child-owner failure
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

Reactive EB source and finite-volume hydro operators now both have direct
owner-only MPI transactions. The root remains one level-wide physics entity
because weighted StateRedist owns overlapping level neighborhoods; its storage
is still synchronized by root tiles. Owner-only molecular transport, sparse
fields, decomposed root-level StateRedist, point-to-point coarse/fine traffic,
topology migration, parallel checkpoint I/O, and a public distributed time
loop remain separate work.

## Owner-only MPI reactive EB AMR molecular transport (`0.113.0`)

- [x] collective transport-record, boundary, switch, interval, and StateRedist consensus
- [x] owner-authoritative start-state synchronization before SSPRK2 blending
- [x] one complete root EB transport Euler update per stage on one physics owner
- [x] ratio-subcycled fine transport on each exclusive child owner
- [x] owner-local coarse/fine time-integrated diffusive flux registers
- [x] deterministic child reflux and root-owner EB-cut conservation closure
- [x] owner-side SSPRK2 blend, EOS recovery, and final reactive average-down
- [x] exact local and global committed Euler-advance accounting
- [x] serial multipatch EB transport state, temperature, and limiter parity
- [x] exact all-rank rollback after a late child-owner failure
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

Reactive EB chemistry, finite-volume hydro, and qualified mixture molecular
transport now all have direct owner-only MPI transactions. Complete EB fields
and topology remain replicated. Sparse rank-local storage, decomposed
root-level StateRedist, point-to-point coarse/fine and halo traffic, dynamic
topology migration, parallel checkpoint I/O, and a public distributed EB AMR
time loop remain separate work.

## Owner-only MPI reactive EB AMR full physics (`0.114.0`)

- [x] outer owner-only `R-T-H-T-R` transaction
- [x] two owner chemistry half-steps
- [x] two owner SSPRK2 transport half-steps
- [x] one owner hydro interval between transport stages
- [x] deferred root, child, limiter, and counter publication
- [x] exact chemistry, hydro, and transport owner accounting
- [x] serial multipatch full-physics state and temperature parity
- [x] exact outer rollback after chemistry and transport precede hydro failure
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The replicated MPI EB AMR bridge now provides individual owner operators and
their complete reactive split composition. Sparse EB fields, point-to-point
halo and coarse/fine traffic, a distributed public timestep loop, dynamic
topology migration, and parallel checkpoint/output remain separate work.

## Sparse MPI reactive EB AMR owner storage (`0.115.0`)

- [x] rank-local sparse root-tile payloads
- [x] rank-local sparse child-patch payloads
- [x] zero numerical allocation for every nonowned entity
- [x] exact local and communicator-wide stored-value accounting
- [x] owner-authoritative scatter from stale replicated inputs
- [x] explicit replicated materialization boundary
- [x] bitwise materialization parity with the established owner synchronizer
- [x] collective invalid-local-payload rejection and output rollback
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

Persistent EB state can now be represented without nonowner numerical
replicas. Owner-only physics still accepts complete temporary fields, so direct
sparse chemistry/hydro/transport, point-to-point materialization, public time
advancement, regridding, checkpointing, and output remain separate work.

## Direct chemistry on sparse MPI reactive EB AMR storage (`0.116.0`)

- [x] direct active-cell chemistry on owned sparse root tiles
- [x] direct active-cell chemistry on owned sparse children
- [x] no nonowner state allocation during reactor execution
- [x] collective interval, tolerance, species, and reaction-width consensus
- [x] exact local and global owner reactor accounting
- [x] post-reaction materialize, root-owner average-down, and sparse re-scatter
- [x] bitwise serial patch-set chemistry parity
- [x] exact sparse rollback after a late child-owner reactor failure
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

Chemistry is the first 2D EB operator to consume persistent sparse owner
payloads directly. Average-down still uses a complete temporary hierarchy;
direct child-to-root synchronization, sparse hydro and transport, public time
advancement, regridding, checkpointing, and output remain separate work.

## Direct average-down on sparse MPI reactive EB AMR storage (`0.117.0`)

- [x] child-owner volume-fraction-weighted conserved restriction
- [x] one coarse-footprint buffer broadcast per child
- [x] root-tile-owner application over exact child intersections
- [x] covered coarse-cell preservation
- [x] owner-local EOS temperature recovery
- [x] no complete temporary root or child hierarchy
- [x] unchanged sparse input/output allocation count
- [x] bitwise parity with replicated reactive average-down
- [x] exact rollback after a late child restriction fails EOS recovery
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

Sparse chemistry now remains below the complete materialization boundary from
input through commit. Restriction buffers are still communicator broadcasts;
point-to-point child/root traffic, sparse hydro and transport, public time
advancement, regridding, checkpointing, and output remain separate work.

## Sparse MPI reactive EB AMR full-physics transaction (`0.118.0`)

- [x] sparse input and output for the complete `R-T-H-T-R` split
- [x] first direct sparse chemistry half-step
- [x] one complete temporary hierarchy across the central `T-H-T` window
- [x] immediate owner re-scatter before final direct sparse chemistry
- [x] unchanged committed sparse allocation count
- [x] deferred chemistry, hydro, transport, and limiter publication
- [x] exact owner and global operator accounting
- [x] serial multipatch full-physics state, temperature, and limiter parity
- [x] exact outer rollback when hydro rejects after chemistry and transport
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The public sparse transaction now covers the complete qualified physics split.
Hydro and transport still operate inside one replicated compatibility window;
direct sparse hydro/transport, point-to-point traffic, time-loop control,
regridding, checkpointing, and output remain separate work.

## Direct hydro on sparse MPI reactive EB AMR storage (`0.119.0`)

- [x] root-tile assembly into a level-wide StateRedist temporary
- [x] one root-owner hydro update and root flux synchronization
- [x] owner-local child exterior construction and ratio subcycling
- [x] owner-local coarse/fine flux-register accumulation and reflux
- [x] no allocation or synchronization of nonowner child payloads
- [x] direct corrected-root scatter and sparse reactive average-down
- [x] unchanged committed sparse allocation count
- [x] exact local and global level-advance accounting
- [x] serial multipatch hydro state and temperature parity
- [x] exact rollback after a late final-child failure
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

Fine-level hydro now remains owner-local throughout the transaction. Root state
and flux arrays are still level-wide broadcasts because root EB StateRedist is
not decomposed; root halo/redistribution decomposition, targeted traffic,
time-loop control, regridding, checkpointing, and output remain separate work.

## Direct transport on sparse MPI reactive EB AMR storage (`0.120.0`)

- [x] SSPRK2 transport transaction over sparse persistent state
- [x] root-tile assembly without fine-child materialization
- [x] root-owner transport, StateRedist, and flux synchronization per stage
- [x] owner-local child exterior construction and ratio subcycling
- [x] owner-local diffusive flux-register accumulation and reflux
- [x] distributed composite integral for EB-cut conservation closure
- [x] no allocation or synchronization of nonowner child payloads
- [x] direct sparse average-down after every Euler stage and final blend
- [x] unchanged committed sparse allocation count
- [x] exact local and global Euler-advance accounting
- [x] serial state, temperature, and limiter parity
- [x] exact rollback after a late final-child failure
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

Fine-level transport now remains owner-local through both SSPRK2 stages. Root
state, temperature, and flux arrays are still level-wide broadcasts because
root transport and EB StateRedist are not decomposed. The outer sparse
`R-T-H-T-R` driver still uses its established central compatibility window;
wiring these direct component operators into that transaction, targeted
root/coarse-fine traffic, time-loop control, regridding, checkpointing, and
output remain separate work.

## End-to-end full physics on sparse MPI reactive EB AMR storage (`0.121.0`)

- [x] direct sparse `R-T-H-T-R` composition
- [x] direct sparse chemistry for both reaction half-steps
- [x] direct sparse SSPRK2 transport for both transport half-steps
- [x] direct sparse hydro for the central hyperbolic interval
- [x] no complete fine-child hierarchy between physics operators
- [x] unchanged committed sparse allocation count
- [x] exact local and global chemistry, hydro, and transport accounting
- [x] serial root and child state and temperature parity
- [x] serial transport-limiter parity
- [x] outer rollback after first chemistry and transport succeed
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The complete qualified physics split now retains globally single-copy fine
payloads from sparse input through sparse output. Component entrypoints still
assemble and synchronize level-wide root temporaries. Root decomposition,
targeted root and coarse/fine traffic, public time-loop control, regridding,
checkpointing, and output remain separate work.

## Targeted sparse EB average-down traffic (`0.122.0`)

- [x] one child-owner restriction buffer per fine child
- [x] exact intersecting root-owner recipient discovery
- [x] duplicate recipient elimination when one owner has multiple root tiles
- [x] no restriction allocation or receive on unrelated ranks
- [x] point-to-point transfer only to distinct remote recipients
- [x] local application when child and root tile share an owner
- [x] unchanged covered-cell preservation and owner-local EOS recovery
- [x] exact local and communicator transfer accounting
- [x] bitwise serial chemistry and average-down parity
- [x] exact state and published-count rollback after EOS rejection
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

Sparse average-down no longer broadcasts every child restriction across the
communicator. Root physics fields, coarse-time exterior state, fluxes, and
reflux corrections still use level-wide synchronization in the direct hydro
and transport operators. Root decomposition, targeted physics traffic, public
time-loop control, regridding, checkpointing, and output remain separate work.

## Targeted direct sparse hydro root traffic (`0.123.0`)

- [x] one packed root-tile gather per non-root physics owner tile
- [x] root physics state and flux allocation only on the root owner
- [x] one packed root bundle per distinct remote child owner
- [x] no full root allocation on unrelated ranks
- [x] one correction payload in each direction per remote child reflux
- [x] one packed final row-band scatter per remote root tile
- [x] no all-rank hydro numerical field broadcast
- [x] exact local and communicator payload-transfer accounting
- [x] unchanged sparse allocation count after commit
- [x] serial root and child state and temperature parity
- [x] exact state, advance-count, and transfer-count rollback
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

Hydro still advances the complete root level on one physics owner because root
EB reconstruction and StateRedist are not decomposed, but its full numerical
arrays no longer reach unrelated ranks. Root decomposition, targeted
transport traffic, public time-loop control, regridding, checkpointing, and
output remain separate work.

## Targeted direct sparse transport root traffic (`0.124.0`)

- [x] one packed root-tile gather per Euler stage and non-root owner tile
- [x] root transport state, temperature, RHS, and flux allocation only on the
  root physics owner and actual child owners
- [x] one packed root bundle per Euler stage and distinct remote child owner
- [x] one correction payload in each direction per stage and remote child
- [x] one packed corrected row-band scatter per Euler stage and remote tile
- [x] targeted two-candidate gather and final SSPRK2 row-band scatter
- [x] tile-local EB-cut conservation closure with only an `nvar` vector
  broadcast
- [x] no all-rank transport root-field broadcast or unrelated full-root
  allocation
- [x] exact local and communicator payload-transfer accounting
- [x] unchanged sparse allocation count after commit
- [x] serial root, child, temperature, and limiter parity
- [x] exact state, Euler-count, limiter, and transfer-count rollback
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

Transport still advances the complete root level on one physics owner because
the root diffusive stencil and StateRedist remain level-wide. Full root fields
now exist only on that owner and ranks that actually own fine children. Root
physics decomposition, public time-loop control, regridding, checkpointing,
and output remain separate work.

## Sparse owner-local EB timestep selection (`0.125.0`)

- [x] one targeted sparse root gather to the root physics owner
- [x] root EB hyperbolic CFL evaluation only on that owner
- [x] optional root parabolic transport limit only on that owner
- [x] hydro and transport limits evaluated only for locally owned children
- [x] refinement-ratio scaling from fine dt to coarse interval
- [x] one communicator-minimum global stable timestep
- [x] exact root-gather send accounting
- [x] serial patch-set hydro/transport timestep parity
- [x] finite negative child-state collective rejection
- [x] unchanged sparse state and zero dt/transfer publication on rejection
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The distributed sparse hierarchy can now choose its own qualified coarse
interval without replicated fine payloads or all-rank root fields. A public
multi-step driver, dynamic topology, checkpointing, and output remain separate
work.

## Public sparse MPI reactive EB AMR time loop (`0.126.0`)

- [x] public target-time loop over the direct sparse `R-T-H-T-R` transaction
- [x] fresh owner-local hydro/transport stability selection before every step
- [x] exact final-time clipping with a positive finite accepted interval
- [x] communicator consensus for clock, limits, tolerances, and controls
- [x] time and total-step publication only after a complete split-step commit
- [x] committed-only chemistry, hydro, transport, limiter, and timestep-traffic
  accounting
- [x] preservation of earlier accepted states and diagnostics at a later
  total-step limit
- [x] serial dynamic-timestep root, child, and temperature parity
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The static sparse EB hierarchy can now advance its own full-physics clock
without replicated fine payloads. Dynamic sparse topology, checkpoint/restart,
parallel output, and decomposition of the level-wide root physics kernel
remain separate work.

## Transactional sparse MPI reactive EB AMR topology rebuild (`0.127.0`)

- [x] public explicit-plan sparse child-topology replacement
- [x] serial EB average-down, PCM initialization, and old/new overlap retention
- [x] deterministic subcycle-weighted owner recomputation for the new topology
- [x] return to rank-local root tiles and exclusive child payloads after regrid
- [x] exact one-copy global stored-value accounting after owner redistribution
- [x] atomic distribution, sparse state, and geometry-template publication
- [x] invalid refinement-ratio rejection before materialization or mutation
- [x] exact root, child, temperature, and topology parity with serial regrid
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The explicit regrid removes the fixed-topology API restriction but currently
uses a replicated compatibility window during topology changes. Tag planning,
scheduled clock integration, direct old/new overlap transfer,
checkpoint/restart, parallel output, and root physics decomposition remain
separate work at this milestone.

## Scheduled tag-driven sparse MPI reactive EB AMR regrid (`0.128.0`)

- [x] root-owner-only temperature-gradient tag evaluation and collection plan
- [x] compact ordered-plan metadata broadcast without root numerical fields
- [x] public caller geometry-builder callback for planned EB child rectangles
- [x] accepted-step cadence integrated into the sparse target-time loop
- [x] atomic physics-step plus distribution/state/template regrid commit
- [x] committed-only regrid-evaluation, topology-change, and traffic accounting
- [x] independent serial dynamic-timestep and scheduled-regrid reference
- [x] exact root, child, temperature, topology, and limiter parity
- [x] geometry-builder failure rollback of state, clock, hierarchy, and counts
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The sparse EB clock now owns a complete two-level dynamic-topology lifecycle.
At this milestone regrid events still use the `0.127.0` replicated
compatibility window after the owner-local plan is broadcast. Direct overlap
migration, sparse checkpoint/restart and parallel output, arbitrary-depth
dynamic EB topology, and decomposition of the level-wide root physics kernel
remain separate work.

## Direct sparse MPI reactive EB AMR topology transfer (`0.129.0`)

- [x] old fine-to-root average-down through targeted restriction buffers
- [x] averaged root-tile assembly only on distinct new child owners
- [x] owner-local PCM initialization of every new child
- [x] local or point-to-point same-ratio old/new overlap rectangle retention
- [x] owner-local active-cell temperature recovery after overlap migration
- [x] no all-rank root or child numerical-field materialization during regrid
- [x] exact restriction, prolongation, and overlap send accounting
- [x] transfer counters exposed through explicit, tagged, and scheduled APIs
- [x] candidate-only mutation and atomic distribution/state/template commit
- [x] late overlap-geometry mismatch rollback after restriction and PCM staging
- [x] exact serial root, child, temperature, and topology parity
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The two-level sparse EB lifecycle now stays globally single-copy through both
physics and dynamic topology changes. Geometry and compact topology metadata
remain replicated intentionally. Sparse checkpoint/restart and parallel output,
arbitrary-depth dynamic EB topology, and decomposition of the level-wide root
physics kernel remain separate work.

## Root-only sparse MPI reactive EB AMR materialization (`0.130.0`)

- [x] public caller-selected root gather for sparse root tiles and fine children
- [x] one packed state/temperature send per remotely owned entity
- [x] no complete root or patch-set field allocation on non-root ranks
- [x] exact root and every-child parity with established all-rank materialization
- [x] independent local-sender and communicator-summed transfer accounting
- [x] rank-consistent root selection and collective sparse-input validation
- [x] unallocated outputs and zero published traffic on invalid input
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

Checkpoint and output adapters can now consume a complete hierarchy on one
writer rank without replicating numerical fields on every rank. The existing
formatted checkpoint and CSV writers are not yet connected to this boundary;
rank-independent restart distribution and parallel file formats remain
separate work.

## Sparse MPI reactive EB AMR checkpoint and CSV writers (`0.131.0`)

- [x] public root-selected sparse formatted checkpoint writer
- [x] public root-selected sparse root/child CSV writer
- [x] complete numerical arrays allocated only inside the selected writer
- [x] existing serial multipatch checkpoint schema retained without conversion
- [x] checkpoint read-back parity for root and every child field
- [x] one nonempty root CSV and one nonempty deterministic CSV per child
- [x] collective writer status after root-local filesystem operations
- [x] exact successful gather traffic and zero published traffic on I/O failure
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The adapter does not yet read and redistribute a checkpoint directly into a
new sparse owner map. Parallel file formats and atomic multi-file CSV rename
remain separate work.

## Root-only sparse MPI reactive EB AMR checkpoint restart (`0.132.0`)

- [x] public root-to-owner direct sparse hierarchy scatter
- [x] public root-only formatted checkpoint read adapter
- [x] no complete numerical read allocation or broadcast on non-root ranks
- [x] one packed send per remote root tile or fine child
- [x] exact root-local and communicator-summed restart transfer accounting
- [x] root and child owner-local field parity after checkpoint round trip
- [x] collective clock, step, regrid, minimum-timestep, and density metadata
- [x] empty sparse state, zero metadata, and zero traffic on read failure
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

The caller still supplies a replicated geometry and full patch-set template
whose topology must match the checkpoint. A geometry-only topology descriptor,
checkpoint-driven dynamic topology reconstruction, and a cross-run rank-count
redistribution application remain separate work.

## Geometry-only sparse MPI reactive EB AMR restart topology (`0.133.0`)

- [x] public child geometry/topology descriptor without state or temperature
- [x] direct construction from fine geometries and an ordered patch collection
- [x] validated extraction from an established full reactive patch set
- [x] distribution and sparse-payload validation against geometry only
- [x] direct root-to-owner restart scatter without non-root full child fields
- [x] root-only checkpoint read accepting the geometry-only descriptor
- [x] compatibility wrappers retaining the former full patch-set interfaces
- [x] exact owner-local root/child parity and restart transfer accounting
- [x] collective empty-state and zero-traffic rollback for invalid topology
- [x] OpenMPI one-, two-, four-, and eight-rank gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

EB geometry arrays and compact patch boxes remain intentionally replicated,
but restart no longer needs replicated child conserved or temperature fields.
Checkpoint-driven topology reconstruction, geometry-only conversion of the
writer and physics/regrid compatibility boundaries, arbitrary-depth dynamic EB
topology, and decomposition of the level-wide root physics kernel remain
separate work.

## Arbitrary-depth geometry-only reactive EB patch-tree topology (`0.134.0`)

- [x] root plus a runtime-sized sequence of refinement relations
- [x] multiple ordered parents and children on every refinement level
- [x] parent/local-child to flattened child index mapping
- [x] per-relation refinement ratios and exact EB patch reconstruction
- [x] parent geometry, nesting, and separated-sibling validation
- [x] direct four-level, two-branch construction gate
- [x] whole-tree transactional rebuild and exact no-op detection
- [x] invalid-parent rollback preserving the accepted topology
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

This milestone removes the fixed two-/three-level limit from geometry topology
metadata only. Arbitrary-depth numerical fields, conservative state migration,
owner distribution, physics recursion, timestep selection, checkpoint I/O, and
dynamic tagging remain separate work.

## Arbitrary-depth reactive EB patch-tree state migration (`0.135.0`)

- [x] numerical state and temperature node for every runtime topology patch
- [x] parent-first PCM initialization across arbitrary depth and branching
- [x] deepest-first reactive average-down across every tree relation
- [x] complete composite conserved-vector integral at arbitrary depth
- [x] transactional whole-tree numerical rebuild from a collapsed old root
- [x] same-resolution physical overlap retention with local EB metric checks
- [x] parent-first regeneration after retained overlap for deeper descendants
- [x] active-cell NASA7 temperature recovery before candidate publication
- [x] exact no-op behavior for an identical topology plan
- [x] invalid-plan rollback preserving topology, state, and temperature
- [x] four-level, two-branch and shifted-patch conservation/migration gates
- [x] GNU Fortran Release and bounds/FPE-checked Debug qualification

This milestone removes the fixed depth limit from serial numerical hierarchy
storage, synchronization, and conservative topology migration. Runtime physics
recursion, owner distribution, distributed migration, timestep selection,
checkpoint I/O, and dynamic tagging remain separate work.

## MPI owner-tiled reactive EB root hydro (`0.136.0`)

- [x] one bounded root geometry/state band per local y-tile
- [x] six-row guard covering reconstruction and StateRedist dependencies
- [x] established reactive EB level kernel reused on every tile band
- [x] owned state, temperature, and x-face row publication
- [x] unique shared y-face ownership and complete physical-boundary coverage
- [x] collective root state, temperature, and flux assembly
- [x] no selected-rank full-root hydro advance or result broadcasts
- [x] local tile-advance and bounded computed-cell accounting
- [x] serial root/child parity at one, two, four, and eight ranks
- [x] owner-tiled root hydro retained inside the complete physics split
- [x] collective rollback preserving zero published work after late failure
- [x] OpenMPI Release and bounds/FPE-checked Debug qualification

The replicated owner path still stores the complete root input on every rank
and assembles complete root results for the established child/reflux boundary.
The sparse path removes that input replication in the next milestone while
retaining a root-owner compatibility bundle for child exterior and reflux.

## Sparse MPI owner-tiled reactive EB root hydro (`0.137.0`)

- [x] direct point-to-point state and temperature halo fragments
- [x] no complete sparse root input assembled before hydrodynamic work
- [x] one bounded EB band advance on each root tile owner
- [x] six-row guard shared with the qualified replicated decomposition
- [x] owned input, state, temperature, x-flux, and unique y-face result routing
- [x] complete temporary root bundle only on the established root owner
- [x] unchanged child subcycling, exterior construction, reflux, and scatter
- [x] exact halo, result, child, correction, and scatter transfer accounting
- [x] exact local tile-advance and computed-band-cell accounting
- [x] serial root/child parity at one, two, four, and eight ranks
- [x] sparse full-physics and public time-loop owner accounting
- [x] late-failure rollback with zero published work and traffic

The sparse input remains exclusively tile-owned throughout root hydro. A full
temporary root start/result/flux bundle is still assembled on the root owner
after tile-local computation because fine-child exterior interpolation and
deterministic reflux currently consume level-wide arrays. Sparse root
transport and distributed coarse/fine interface data are subsequent
decomposition work.

## Sparse MPI owner-local EB timestep selection (`0.138.0`)

- [x] exact EB geometry band extraction for every locally owned root tile
- [x] root hyperbolic CFL evaluation directly on exclusive tile state
- [x] root parabolic transport limit directly on exclusive tile state
- [x] fully covered root-tile and child skipping
- [x] child-owner limits retained with refinement-ratio coarse scaling
- [x] one communicator-minimum selection with zero root-field transfers
- [x] exact serial hydro/transport timestep parity
- [x] zero timestep traffic through the public multi-step clock
- [x] scheduled-regrid traffic kept distinct from timestep traffic
- [x] invalid owner state rejection with zero dt and transfer publication

Stable-step selection no longer constructs any complete temporary root field.
The sparse SSPRK2 transport Euler stages still gather root fields to the
selected physics owner and remain the next root decomposition target.

## Sparse MPI owner-local SSPRK2 root blend (`0.139.0`)

- [x] final root conserved-state average directly on exclusive tile owners
- [x] exact EB geometry band extraction for tile-local EOS recovery
- [x] final root temperature recovery directly on exclusive tile owners
- [x] no final-blend start-state gather
- [x] no final-blend second-Euler gather
- [x] no final-blend row-band scatter
- [x] unchanged owner-local child blend and final sparse average-down
- [x] exact four-message root-tile traffic per remote tile across two stages
- [x] serial root/child parity at one, two, four, and eight ranks
- [x] late-failure rollback with zero published work and traffic

The final blend is cell-local and therefore needs no neighboring root values.
The two sparse transport Euler stages still assemble a complete root on the
physics owner for diffusive fluxes, StateRedist, child context, and reflux.
Finite-halo transport stages are the next root decomposition target.

## Sparse MPI owner-tiled SSPRK2 root Euler stages (`0.140.0`)

- [x] direct point-to-point state and temperature halo fragments per stage
- [x] no unconditional selected-root gather before transport work
- [x] one explicit EB transport/StateRedist band per root tile and Euler stage
- [x] six-row guard covering the composed transport and redistribution stencil
- [x] complete-root compatibility band for periodic y-boundary tiles
- [x] owned input, result, temperature, x-flux, and unique y-face routing
- [x] complete temporary root bundle only after distributed tile computation
- [x] unchanged child subcycling, exterior construction, reflux, and scatter
- [x] exact halo, result, child, correction, and scatter transfer accounting
- [x] exact local tile-advance and computed-band-cell accounting
- [x] serial root/child/limiter parity at one, two, four, and eight ranks
- [x] sparse full-physics and public time-loop owner accounting
- [x] late-failure rollback with zero published work and traffic

The sparse root input remains persistently tile-owned, while each target owner
materializes only its required temporary band. Periodic y-boundary targets use
a complete band so cyclic wrap remains exact. Fine-child context and cumulative
reflux also retain a complete temporary root bundle on the established physics
owner. Cyclic band geometry and distributed coarse/fine interface data are
subsequent decomposition work.

## Sparse MPI periodic-edge cyclic transport bands (`0.141.0`)

- [x] boundary-anchored lower- and upper-edge cyclic source-row maps
- [x] physical periodic y boundary retained at the temporary band ends
- [x] six-row transport/StateRedist dependency guard
- [x] one additional row isolating each target guard from the internal gap
- [x] point-to-point assembly split by source owner and contiguous fragment
- [x] absolute EB boundary-centroid y remapping into compact coordinates
- [x] complete-root fallback when the protected footprint spans the root
- [x] exact fragment-transfer and computed-band-cell accounting
- [x] dedicated 14-by-21 root-only cyclic test at one, two, four, and eight ranks
- [x] finite periodic-edge bands exercised at four and eight ranks
- [x] serial state, temperature, and redistribution-limiter parity
- [x] existing dynamic-regrid and public time-loop regression coverage

The cyclic band keeps the serial kernel's periodic boundary at its physical
outer ends instead of placing the wrap at an artificial internal seam. The
post-compute complete root bundle is still assembled for fine-child exterior
data, flux registers, and reflux; distributing that interface remains later
work.

## Compact EB flux-register and sparse correction support (`0.142.0`)

- [x] globally indexed patch-plus-one-cell flux-register storage
- [x] boundary-clipped correction bounds validated with the register
- [x] reflux iteration limited to the stored correction support
- [x] unchanged cut-cell cardinal/diagonal redistribution
- [x] patch-plus-two-cell sparse correction transfer footprint
- [x] cumulative corrections merged in deterministic child order
- [x] unchanged message-count accounting and transactional rollback
- [x] serial flux cancellation, conservation, temperature, and failure gates
- [x] sparse root/child parity at one, two, four, and eight ranks
- [x] public time-loop and scheduled-regrid regression coverage

Each distinct child owner still receives one complete root start/end/flux
input bundle per Euler stage. Compact exterior and interface-flux routing are
the next boundary; this milestone removes only the repeated full-root reflux
correction round trips.

## Compact sparse child transport context (`0.143.0`)

- [x] raw start/end coarse exterior samples stored only on four child edges
- [x] compact-context reconstruction identical to the complete-root builder
- [x] coarse interface flux accumulated into the compact register on root
- [x] child owner advances only fine state and compact register data
- [x] evolved fine state and accumulated mismatch returned in one payload
- [x] deterministic reflux performed on the root physics owner
- [x] corrected fine state returned without any complete root array
- [x] exact three-message remote-child transaction per Euler stage
- [x] context payload required to be smaller than the former root bundle
- [x] serial root/child parity retained at one, two, four, and eight ranks

The complete root result still exists transiently on the root physics owner
after owner-tiled Euler work. Distributing final reflux across root-tile owners
is subsequent work; transport child owners no longer receive that result.

## Compact child-local reactive reflux support (`0.144.0`)

- [x] nonreactive reflux accepts globally indexed coarse support arrays
- [x] reactive reflux recovers only support and fine active-cell temperatures
- [x] complete-root reflux remains an exact support-wrapper compatibility API
- [x] patch-plus-two support routed with compact child exterior context
- [x] fine child state remains exclusively on its owner through reflux
- [x] only corrected coarse support returns to the root physics owner
- [x] deterministic overlapping-support merge remains in child order
- [x] two messages per remote child per transport Euler stage
- [x] combined context/support payload smaller than the former root bundle
- [x] strictly smaller globally indexed support/full reactive reflux unit parity

The root physics owner still assembles the temporary complete root Euler
result and merges returned supports before row scatter. Direct context and
correction routing between root-tile owners is subsequent work.

## Compact coarse interface-flux accumulation (`0.145.0`)

- [x] coarse flux accumulation accepts globally indexed x/y face rectangles
- [x] validation requires every active coarse/fine interface face
- [x] invalid or nonfinite support leaves the register bitwise unchanged
- [x] complete-root accumulation remains a thin compatibility wrapper
- [x] compact/full correction arrays are bitwise identical
- [x] sparse MPI transport passes only each child's interface rectangles

The root physics owner still owns the complete temporary root x/y flux bundle
after tile computation. This API removes that shape requirement from flux-
register accumulation so later root-tile-to-child routing can supply only the
interface fragments.

## Direct root-tile coarse-flux routing (`0.146.0`)

- [x] each root tile retains its computed x-flux rows and unique y-faces
- [x] only child-intersecting face fragments route to the child owner
- [x] child-owner coverage checks reject missing ownership gaps
- [x] compact coarse register is accumulated on the child owner
- [x] root-to-child state context contains no flux-register payload
- [x] context plus interface-flux values remain below the legacy root bundle
- [x] exact direct-fragment message accounting for both SSPRK2 Euler stages

The root physics owner still assembles the temporary complete root result and
flux bundle for exterior extraction, support state, boundary closure, ordered
support merge, and final row scatter. Direct state/support routing is the next
ownership boundary.

## Compact exterior state-context support (`0.147.0`)

- [x] context extraction accepts globally indexed coarse support arrays
- [x] patch-plus-one start/end state and temperature is sufficient
- [x] complete-root extraction remains a support-wrapper compatibility API
- [x] compact/full context reconstruction is bitwise identical
- [x] incomplete, out-of-root, or nonfinite support rejects transactionally

Sparse MPI still extracts state context on the root physics owner in this
milestone. The support API removes the complete-root shape dependency required
before root tile owners can route start/end state fragments directly.

## Direct root-tile state/support routing (`0.148.0`)

- [x] tile owners retain start, uncorrected-end, and current corrected state
- [x] patch-plus-two state/temperature fragments route directly to child owner
- [x] child owner extracts the established four-edge exterior context locally
- [x] cumulative corrected support remains visible in deterministic child order
- [x] corrected fragments return directly to intersecting root tile owners
- [x] final corrected root rows commit without a root-owner scatter
- [x] state plus interface-flux payload remains below the legacy root bundle
- [x] exact three-route tile/child message accounting covers both Euler stages

The root physics owner still receives the complete temporary tile result and
flux bundle for compatibility validation and cut-boundary flux closure. This
milestone removes it from child state/reflux routing but does not yet remove
that complete post-compute assembly or claim a measured speedup.

## Owner-local root transport result (`0.149.0`)

- [x] transport Euler publishes only owner-local tile state and flux records
- [x] no remote tile sends a complete owned result to the root physics owner
- [x] no complete root transport state, temperature, or flux array is allocated
- [x] left/right boundary flux change is accumulated on every tile owner
- [x] lower/upper boundary flux change is owned only by the edge tiles
- [x] one communicator sum replaces root-array boundary closure and broadcast
- [x] exact transport traffic contains only halos and direct child fragments
- [x] root-only and cut-interface paths retain serial numerical parity gates

Complete hierarchy materialization remains available only at explicit
checkpoint, output, restart, and legacy compatibility boundaries. Hydro is
unchanged and still assembles its temporary complete-root physics bundle.

## Compact sparse hydro child context (`0.150.0`)

- [x] four-edge hydro start/end context extracted on the root physics owner
- [x] current patch-plus-two corrected state/temperature sent per child
- [x] only child-intersecting coarse x/y flux rectangles are transferred
- [x] compact state and flux payload packed into one remote-child message
- [x] coarse flux register accumulated on the child owner from global bounds
- [x] hydro reflux runs on patch-plus-two support on the child owner
- [x] only corrected support returns to the root owner
- [x] remote child owners allocate no complete root hydro or flux field
- [x] deterministic child ordering and final root scatter remain unchanged

The root physics owner still assembles the owner-tiled complete hydro result
and scatters final corrected rows. Direct tile-to-child state/flux routing and
tile-owner correction are later ownership boundaries.

## Direct hydro coarse-flux routing (`0.151.0`)

- [x] hydro tile owners retain x rows and uniquely owned y-faces
- [x] tile-to-root hydro result payload omits all flux values
- [x] root-to-child state context contains no flux values
- [x] intersecting tile owners route compact flux fragments to each child
- [x] child owner verifies complete x/y interface coverage before accumulation
- [x] root physics owner allocates no complete hydro x/y flux field
- [x] exact hydro traffic counts one message per remote tile/child intersection
- [x] complete root state assembly and corrected-row scatter remain explicit

The next ownership boundary is hydro state/support routing directly between
root tile owners and child owners, including ordered correction return. This
milestone makes no measured speedup claim.

## Owner-local hydro result and direct state/support routing (`0.152.0`)

- [x] hydro tile owners retain start, end, and current corrected state
- [x] patch-plus-two state fragments route directly to each child owner
- [x] child owner extracts the four-edge hydro context from assembled support
- [x] ordered reflux corrections return directly to intersecting tile owners
- [x] final corrected hydro root rows commit locally
- [x] no complete root hydro state, temperature, or flux result is allocated
- [x] no tile-result or final corrected-row scatter message remains
- [x] exact hydro traffic contains only halos and direct child fragments
- [x] superseded private root-bundle communication helpers are removed

Sparse hydro and transport now share the same owner-local state, flux, and
correction routing boundary. Complete hierarchy materialization remains only
for explicit output, checkpoint, restart, and legacy compatibility APIs. This
milestone makes no measured speedup claim.

## Arbitrary-depth reactive EB patch-tree timestep (`0.153.0`)

- [x] one shared active-cell EB CFL kernel below driver and AMR layers
- [x] every root and child patch evaluated on its own EB geometry
- [x] fully covered nodes skipped while an entirely inactive tree rejects
- [x] cumulative refinement-product conversion to one root interval
- [x] arbitrary-depth and branching traversal with no fixed level count
- [x] species-layout, finite-CFL, node-conversion, and scale validation
- [x] deterministic zero timestep on every rejected path
- [x] read-only hierarchy contract on accepted and rejected paths
- [x] four-level two-branch reference reduction gate
- [x] forced deepest-node limiting gate with exact subcycle scaling
- [x] nonfinite-CFL rollback gate

This milestone removes timestep selection from the remaining serial patch-tree
gaps. Arbitrary-depth hydrodynamics, chemistry, molecular transport, public
clock ownership, dynamic tagging, checkpoint I/O, and MPI ownership for this
new 2D EB tree remain separate work.

## Arbitrary-depth reactive EB patch-tree hydrodynamics (`0.154.0`)

- [x] runtime-depth recursive EB level advancement
- [x] multiple ordered children attached to their actual parent patch
- [x] exact per-relation ratio subcycling and parent-time interpolation
- [x] one independent coarse/fine EB flux register per child
- [x] deterministic reflux and reactive average-down at every relation
- [x] composite subtree conservation closure against outer-boundary flux
- [x] active unrefined-parent correction with NASA7 temperature recovery
- [x] deepest-first final synchronization before atomic tree commit
- [x] optional committed-only per-level node-advance counts
- [x] four-level two-branch schedule `[1, 4, 8, 8]`
- [x] mass, total-energy, and every-species composite conservation gates
- [x] fixed three-level field/temperature parity and schedule `[1, 2, 4]`
- [x] invalid-solver rollback preserving every node and zeroing counts

This milestone removes hyperbolic advancement from the remaining serial
patch-tree gaps. Arbitrary-depth chemistry, molecular transport, full-physics
composition, public clock ownership, dynamic tagging, checkpoint I/O, and MPI
ownership for this new 2D EB tree remain separate work.

## Arbitrary-depth reactive EB patch-tree chemistry (`0.155.0`)

- [x] EB active-mask chemistry on every runtime tree patch
- [x] standalone transactional chemistry with deepest-first synchronization
- [x] `R-H-R` Strang composition around recursive tree hydrodynamics
- [x] one private candidate across both chemistry halves and hydro
- [x] committed-only chemistry and hydro per-level call counts
- [x] four-level two-branch chemistry schedule `[1, 2, 2, 1]`
- [x] four-level two-branch Strang schedules `[2, 4, 4, 2]` and `[1, 4, 8, 8]`
- [x] composite mass, total-energy, species-closure, and activity gates
- [x] fixed three-level Strang field and temperature parity
- [x] exact rollback after chemistry followed by invalid hydrodynamics

This milestone removes active-cell reaction and hydro/chemistry Strang
composition from the remaining serial patch-tree gaps. Arbitrary-depth
molecular transport, `R-T-H-T-R` full physics, public clock ownership, dynamic
tagging, checkpoint I/O, and MPI ownership remain separate work.

## Arbitrary-depth reactive EB patch-tree transport (`0.156.0`)

- [x] runtime-depth recursive EB transport Euler advancement
- [x] parent start/end interpolation and exact per-relation child subcycling
- [x] one independent diffusive EB flux register per ordered child
- [x] deterministic reflux, average-down, and subtree conservation closure
- [x] two complete recursive Euler stages with node-wise SSPRK2 blending
- [x] EOS temperature recovery and deepest-first final synchronization
- [x] committed-only minimum limiter theta and per-level node counts
- [x] fixed three-level SSPRK2 parity with schedule `[2, 4, 8]`
- [x] four-level branching schedule `[2, 4, 16, 16]`
- [x] composite conservation, changed-state, validity, and rollback gates

This milestone removes standalone molecular transport from the remaining
serial patch-tree gaps. Combined `R-T-H-T-R` full physics, public clock
ownership, dynamic tagging, checkpoint I/O, and MPI ownership remain separate
work.

## Arbitrary-depth reactive EB patch-tree full physics (`0.157.0`)

- [x] one private candidate across `R-T-H-T-R` split physics
- [x] qualified all-node chemistry on both reaction half-steps
- [x] qualified recursive SSPRK2 transport on both transport half-steps
- [x] qualified recursive hydrodynamics on the complete interval
- [x] committed-only chemistry, transport, and hydro per-level counts
- [x] committed-only minimum transport limiter theta
- [x] three-level fixed-path field and temperature parity
- [x] three-level schedules `[2, 2, 2]`, `[4, 8, 16]`, and `[1, 2, 4]`
- [x] four-level branching schedules `[2, 2, 4, 2]`, `[4, 8, 32, 32]`, and
  `[1, 2, 8, 8]`
- [x] composite conservation, thermodynamic validity, and late rollback gates

This milestone removes split full-physics composition from the remaining
serial patch-tree gaps. Dynamic tagging, checkpoint I/O, and MPI ownership
remain separate work.

## Public arbitrary-depth reactive EB patch-tree time loop (`0.158.0`)

- [x] combined all-node hydro and explicit transport stable-step selection
- [x] cumulative refinement scaling on every runtime level and branch
- [x] exact stop-time clipping with caller-owned time and total step count
- [x] one private tree candidate per complete `R-T-H-T-R` step
- [x] committed-only minimum timestep and minimum transport limiter theta
- [x] committed-only accumulated chemistry, transport, and hydro level counts
- [x] exact parity with an independently repeated two-step reference sequence
- [x] maximum-step termination retaining the last committed prefix
- [x] exact first-step failure rollback with zero public accounting

This milestone removes public clock ownership from the remaining serial
patch-tree gaps. Dynamic tagging, checkpoint I/O, and MPI ownership remain
separate work.

## MPI arbitrary-depth reactive EB patch-tree ownership (`0.159.0`)

- [x] one deterministic owner for every runtime level/patch pair
- [x] configurable cumulative-subcycle work weighting
- [x] exact per-rank allocated-cell, entity, and weighted-work accounting
- [x] collective topology geometry, relation, and control consensus
- [x] owner-authoritative state and temperature publication
- [x] one private replicated candidate and all-rank commit boundary
- [x] four-level branching coverage at one, two, four, and eight ranks
- [x] collective rejection of rank-local invalid state and control mismatch

This milestone establishes ownership while retaining replicated field
allocation. Sparse owner storage, direct owner migration, distributed
timestep reduction, and owner-local recursive physics remain separate work.

## MPI sparse arbitrary-depth reactive EB patch-tree storage (`0.160.0`)

- [x] replicated topology and owner metadata with owner-only field allocation
- [x] exact local allocated-cell and node accounting
- [x] explicit owner-to-replicated materialization boundary
- [x] direct old-owner to new-owner state and temperature migration
- [x] local copy for unchanged ownership and no nonowner field allocation
- [x] one private sparse migration candidate and all-rank commit boundary
- [x] exact field parity before and after rotated ownership
- [x] collective invalid-map rejection with zero transfer accounting
- [x] one-, two-, four-, and eight-rank Debug and Release coverage

This milestone removes replicated persistent fields from the arbitrary-depth
MPI tree. Distributed sparse timestep reduction and owner-local recursive
physics remain separate work.

## MPI owner-local arbitrary-depth EB patch-tree timestep (`0.161.0`)

- [x] owner-local hyperbolic EB stability evaluation
- [x] owner-local explicit mixture-transport stability evaluation
- [x] cumulative refinement-product conversion to one root interval
- [x] communicator-wide minimum with neutral empty-rank contribution
- [x] exact global active-node evaluation accounting
- [x] collective CFL, transport-flag, distribution, and sparse-state checks
- [x] exact parity with the complete serial patch-tree selector
- [x] collective rank-local control-mismatch rejection with neutral outputs
- [x] one-, two-, four-, and eight-rank Debug and Release coverage

This milestone removes materialization from distributed timestep selection.
Owner-local recursive hydro, transport, chemistry, and public clock execution
remain separate work.

## MPI owner-local arbitrary-depth EB patch-tree chemistry (`0.162.0`)

- [x] owner-local active-cell chemistry on one private sparse candidate
- [x] owner-local conserved-state temperature recovery
- [x] per-node collective acceptance and committed-only advance accounting
- [x] deepest-first parent/child synchronization
- [x] local restriction for shared ownership
- [x] direct child-state transfer to a distinct parent owner
- [x] exact map-derived restriction-transfer accounting
- [x] exact state and temperature parity with the complete serial tree
- [x] collective control-mismatch rollback with zero public accounting
- [x] one-, two-, four-, and eight-rank Debug and Release coverage

This milestone removes materialization from chemistry and average-down.
Owner-local recursive hydro, transport, and public clock execution remain
separate work.

## MPI sparse arbitrary-depth EB composite integrals (`0.163.0`)

- [x] whole-tree owner-local composite conserved integral
- [x] arbitrary selected-subtree owner-local composite integral
- [x] direct-child refined masks excluding covered parent cells
- [x] one recursive topology walk without field materialization
- [x] communicator conserved-vector and contributing-node reductions
- [x] exact global subtree-node accounting
- [x] serial integral parity for the whole tree and all five subtrees
- [x] collective valid-selector disagreement rejection with neutral outputs
- [x] one-, two-, four-, and eight-rank Debug and Release coverage

This milestone supplies the before/after conservation reduction needed by
owner-local recursive hydro. Hydro, transport, and public clock execution
remain separate work.

## MPI owner-local arbitrary-depth EB patch-tree hydro (`0.164.0`)

- [x] owner-local recursive EB node advances on one sparse candidate
- [x] compact parent start/end exterior-context transfer per remote edge
- [x] direct fine-flux return per remote child substep
- [x] parent-owner flux-register accumulation and consumption
- [x] direct child-state reflux round trip and ordered average-down
- [x] shared-owner context, flux, reflux, and restriction without traffic
- [x] owner-local subtree conservation closure with prevalidated reductions
- [x] exact map/schedule-derived grouped-transfer accounting
- [x] serial field and composite-conservation parity
- [x] collective control-mismatch rollback with zero public accounting
- [x] one-, two-, four-, and eight-rank Debug and Release coverage

This milestone removes materialization from recursive hydrodynamics. Owner-
local recursive transport and the sparse public full-physics clock remain
separate work.

## MPI owner-local arbitrary-depth EB patch-tree transport (`0.165.0`)

- [x] owner-local recursive transport Euler stages on one sparse candidate
- [x] compact parent start/end exterior-context transfer per remote edge
- [x] direct fine diffusive-flux return per remote child substep
- [x] parent-owner flux-register accumulation and consumption
- [x] direct child-state reflux round trip and ordered average-down
- [x] owner-local StateRedist, SSPRK2 blend, and EOS temperature recovery
- [x] owner-local subtree conservation closure for both Euler stages
- [x] final deepest-first restriction without materialization
- [x] exact map/schedule-derived grouped-transfer accounting
- [x] serial field, limiter, and composite-conservation parity
- [x] collective boundary/control mismatch rollback with zero accounting
- [x] one-, two-, four-, and eight-rank Debug and Release coverage

This milestone removes materialization from recursive molecular transport.
Owner-local full-physics composition and the sparse public clock remain
separate work.

## MPI owner-local arbitrary-depth EB full physics (`0.166.0`)

- [x] one private sparse `R-T-H-T-R` candidate
- [x] optional owner-local chemistry half-steps
- [x] two owner-local SSPRK2 transport half-steps
- [x] one owner-local recursively subcycled hydro interval
- [x] outer consensus before optional-physics branching
- [x] late-stage failure discards every accepted prefix
- [x] committed-only chemistry, transport, and hydro level counts
- [x] separate exact chemistry, transport, and hydro transfer accounting
- [x] minimum transport limiter across both half-steps
- [x] serial full-tree field, limiter, and conservation parity
- [x] collective control-mismatch rollback with zero public accounting
- [x] one-, two-, four-, and eight-rank Debug and Release coverage

This milestone removes materialization from arbitrary-depth sparse full-physics
composition. The owner-local target-time clock remains separate work.

## MPI owner-local arbitrary-depth EB target-time clock (`0.167.0`)

- [x] communicator-consistent accepted/target time and step controls
- [x] fresh owner-local hydro/transport stability reduction before every step
- [x] exact final-interval clipping and target-time publication
- [x] one private sparse full-physics candidate per attempted step
- [x] committed-prefix semantics after a later step failure
- [x] minimum accepted dt and transport limiter tracking
- [x] cumulative timestep-node and per-level physics accounting
- [x] cumulative chemistry, transport, and hydro transfer accounting
- [x] count-overflow rejection before step commit
- [x] collective clock-control mismatch rollback with neutral outputs
- [x] maximum-step rejection without an uncommitted mutation
- [x] serial clock field, limiter, dt, and conservation parity
- [x] one-, two-, four-, and eight-rank Debug and Release coverage

This milestone removes public-clock ownership from the sparse arbitrary-depth
MPI physics path. Dynamic tagging, checkpoint/restart, and output remain
separate lifecycle work.

## Serial arbitrary-depth EB temperature-tagged rebuild (`0.168.0`)

- [x] deepest-first synchronization before tag evaluation
- [x] per-parent active-cell temperature-gradient tagging
- [x] deterministic disconnected-component clustering at every level
- [x] caller-defined EB geometry construction for every planned child
- [x] recursive planning from synchronized root fields to a level ceiling
- [x] graceful branch termination when a parent is too small to tag
- [x] transactional composition with the overlap-preserving tree rebuild
- [x] root-only creation and maximum-depth refinement
- [x] unchanged-plan no-op and tag-free collapse
- [x] composite-conserved-state preservation across topology changes
- [x] geometry-builder rejection with exact accepted-tree rollback
- [x] GNU Fortran Debug and Release coverage in the complete serial suite

This milestone closes serial dynamic topology planning for the arbitrary-depth
2D EB numerical tree. Owner-local MPI planning/migration, checkpoint/restart,
and composite output remain separate lifecycle work.

## MPI owner-local arbitrary-depth EB temperature-tagged rebuild (`0.169.0`)

- [x] tag evaluation only on each prospective parent owner
- [x] compact integer tag-plan reduction without numerical-tree gathering
- [x] caller-defined EB geometry rebuilt consistently on every rank
- [x] deterministic candidate owner-map recomputation
- [x] direct parent-to-new-child PCM initialization
- [x] direct old-owner to new-owner retained-overlap migration
- [x] deepest-first owner-local average-down after topology change
- [x] composite-conserved-state acceptance before atomic commit
- [x] root-only creation through three levels and changed-plan rebuild
- [x] exact unchanged-plan no-op and tag-free collapse
- [x] rank/control/geometry rejection with exact sparse rollback
- [x] exact topology-derived prolongation and restriction traffic checks
- [x] serial field and temperature parity after every accepted rebuild
- [x] one-, two-, four-, and eight-rank Debug and Release coverage

This milestone removes complete numerical-tree materialization from dynamic
arbitrary-depth MPI topology planning and migration. Checkpoint/restart and
composite output remain separate lifecycle work.

## Serial arbitrary-depth EB patch-tree checkpoint (`0.170.0`)

- [x] distinct versioned formatted schema and terminal marker
- [x] ordered species names and lifecycle metadata
- [x] arbitrary-depth parent/child topology records
- [x] complete root and child EB geometry metric storage
- [x] level-major, patch-major conserved state and temperature storage
- [x] bounded levels, patch counts, and geometry allocation sizes
- [x] private candidate reconstruction before publication
- [x] general-EOS temperature recovery from conserved state
- [x] exact branching topology and roundoff-level field round trip
- [x] maximum-depth rejection with neutral outputs
- [x] species-order rejection with neutral outputs
- [x] invalid lifecycle metadata rejection before file replacement
- [x] GNU Fortran Debug and Release coverage in the complete serial suite

This milestone provides the serial storage format required by later root-only
sparse MPI checkpoint I/O. Rank-neutral MPI redistribution and composite output
remain separate lifecycle work.

## Sparse MPI arbitrary-depth EB checkpoint/restart (`0.171.0`)

- [x] selected-root direct gather from numerical-node owners
- [x] serial self-describing checkpoint write only on the I/O root
- [x] checkpoint read only on the I/O root
- [x] compact arbitrary-depth topology and EB geometry broadcast
- [x] deterministic owner-map reconstruction for the current communicator
- [x] caller-selected hyperbolic or parabolic work exponent on restart
- [x] direct root-to-new-owner numerical field scatter
- [x] no stored rank count or owner map
- [x] exact topology-derived gather and scatter transfer accounting
- [x] field, temperature, and lifecycle metadata restart parity
- [x] collective root/depth/exponent/metadata/species consensus
- [x] neutral outputs after collective or file incompatibility
- [x] one-, two-, four-, and eight-rank Debug and Release coverage

This milestone closes rank-neutral checkpoint/restart for the sparse
arbitrary-depth numerical EB tree. Composite hierarchy output remains separate
lifecycle work.

## Arbitrary-depth EB composite CSV output (`0.172.0`)

- [x] one deterministic CSV for a branching tree of arbitrary depth
- [x] direct-child masks excluding every refined parent cell
- [x] one row per composite leaf cell
- [x] level, patch, index, spacing, coordinate, and time columns
- [x] EB volume, classification, boundary length, and normal columns
- [x] conserved and general-EOS primitive/reactive field columns
- [x] selected-root sparse MPI numerical-node gather
- [x] file access and complete tree allocation only on the writer root
- [x] collective root, time, and ordered-species agreement
- [x] exact topology-derived output transfer accounting
- [x] four-level branching serial leaf-count regression
- [x] one-, two-, four-, and eight-rank Debug and Release coverage

This milestone closes the composite diagnostic-output boundary for the serial
and sparse MPI arbitrary-depth numerical EB tree. A runnable arbitrary-depth
2D EB application lifecycle remains separate integration work.

## Runnable serial arbitrary-depth EB application (`0.173.0`)

- [x] dedicated installed `pelef_reactive_eb_patch_tree_2d` executable
- [x] existing reactive-flow, EB-geometry, and AMR namelist reuse
- [x] validated `patch_tree_maximum_levels` control
- [x] root-only initialization with optional initial recursive tagged rebuild
- [x] all-node hydro/transport timestep selection and exact final-time clipping
- [x] transactional `R-T-H-T-R` advance on the public clock
- [x] scheduled arbitrary-depth topology rebuild after committed steps
- [x] scheduled/final self-describing patch-tree checkpoint calls
- [x] self-describing checkpoint restart entrypoint
- [x] final or checkpoint-stop single composite CSV publication
- [x] four populated dynamic levels in the runnable regression case
- [x] CSV topology, spacing, time, EB class, thermodynamic, and species checks
- [x] GNU Fortran Debug and Release coverage in the complete 210-test suite

This milestone connects the serial arbitrary-depth EB numerical tree to a
public input-driven application without changing the established fixed-depth
driver. Application-level checkpoint/restart round-trip comparison and a
sparse MPI public counterpart remain separate qualification work.

## Public patch-tree application restart parity (`0.174.0`)

- [x] uninterrupted four-level dynamic reference application run
- [x] checkpoint write after one committed step and scheduled regrid
- [x] clean application stop immediately after the scheduled write
- [x] separate-process restart through the public executable and namelists
- [x] restored arbitrary-depth geometry, fields, time, and lifecycle counters
- [x] checkpoint magic, schema, species count, level count, and end marker
- [x] stopped output strictly before the requested final time
- [x] reference and restarted outputs at the exact requested final time
- [x] identity-keyed final topology and complete numeric-column parity
- [x] detailed nested full-physics failure context in the public driver
- [x] GNU Fortran Debug and Release coverage in the complete 214-test suite

This milestone closes application-level checkpoint/restart parity for the
serial arbitrary-depth EB tree. Cross-rank public MPI application restart and
input/checkpoint compatibility fingerprints remain later lifecycle work.

## Public sparse-MPI arbitrary-depth EB application (`0.175.0`)

- [x] installed namelist-driven sparse MPI executable
- [x] configurable depth-weighted MPI ownership exponent
- [x] owner-local recursive initial and periodic tagging/regridding
- [x] owner-local timestep and transactional full-physics clock
- [x] sparse checkpoint write and rank-neutral restart entrypoints
- [x] collective integrals and selected-root composite CSV output
- [x] 1/2/4/8-rank four-level topology and complete-field parity
- [x] GNU Fortran Debug and Release MPI qualification
- [x] all 214 serial regressions retained in MPI and non-MPI builds

Fresh startup still constructs one replicated root field before sparse
conversion; arbitrary-depth child fields are never replicated. Removing the
root startup copy and qualifying public cross-rank application restart remain
later lifecycle work.

## Public sparse-MPI application cross-rank restart (`0.176.0`)

- [x] uninterrupted one-rank four-level reference process
- [x] two-rank checkpoint-stop process after one committed root step
- [x] separate four- and eight-rank restart processes
- [x] checkpoint rank count omitted from authoritative restart state
- [x] ownership map recomputed for each restart communicator
- [x] checkpoint written with uniform node weighting
- [x] reference and restarts run with depth-squared weighting
- [x] direct root-to-owner sparse field scatter after root-only read
- [x] exact four-level identity and column agreement after restart
- [x] complete numeric-field parity with the uninterrupted reference
- [x] GNU Fortran Debug and Release MPI qualification
- [x] all 214 serial regressions retained in MPI and non-MPI builds

This milestone closes the public cross-rank application restart composition.
Fresh initialization still constructs one temporary replicated root field
before converting it to owner-tiled sparse storage; removing that startup copy
is the next public sparse-lifecycle boundary.

## Owner-local public sparse-MPI startup (`0.177.0`)

- [x] root geometry and topology built before numerical initialization
- [x] deterministic root-node ownership selected before field allocation
- [x] reactive root state initialized on exactly one owning rank
- [x] root state and temperature remain unallocated on every non-owner
- [x] collective state-width, shape, finiteness, and temperature validation
- [x] allocatable field ownership moved directly into the sparse root node
- [x] no numerical root broadcast or replicated tree construction
- [x] source state and temperature deallocated by ownership transfer
- [x] existing recursive initial tagging begins from the sparse root
- [x] unchanged four-level 1/2/4/8-rank composite parity
- [x] unchanged cross-rank checkpoint/restart parity
- [x] GNU Fortran Debug and Release MPI qualification
- [x] all 214 serial regressions retained in MPI and non-MPI builds

This milestone removes the last intentionally replicated numerical field from
fresh public sparse-MPI application startup. Geometry and relation metadata
remain replicated for deterministic planning; application/checkpoint input
compatibility fingerprints remain later lifecycle work.

## Public patch-tree checkpoint compatibility fingerprint (`0.178.0`)

- [x] structured mesh, EB, physics, and regridding fingerprint
- [x] public serial and sparse-MPI schema-2 checkpoint writes
- [x] fingerprint validation before topology or field payload reads
- [x] exact integer, flag, and method-name compatibility checks
- [x] round-trip-safe real-control compatibility checks
- [x] final-time, output, checkpoint schedule, rank, and ownership mutability
- [x] serial incompatible-CFL restart rejection with neutral outputs
- [x] collective MPI incompatible-CFL restart rejection
- [x] unchanged 2-to-4/eight-rank valid restart parity
- [x] schema-1 low-level verification compatibility retained
- [x] GNU Fortran Debug and Release MPI qualification
- [x] all 215 serial regressions retained in MPI and non-MPI builds

This milestone closes explicit public application/checkpoint input
compatibility for the arbitrary-depth EB tree.

## Interface-local multilevel EB conservation closure (`0.179.0`)

- [x] topology-derived coarse/fine interface support for every direct child
- [x] clipped three-by-three parent recipient neighborhoods
- [x] exclusion of refined and EB-covered parent cells
- [x] union support for sibling child rectangles
- [x] fixed three-level hydro and transport integration
- [x] serial and sparse-MPI multipatch integration
- [x] serial and sparse-MPI arbitrary-depth tree integration
- [x] fluid-volume normalization and EOS recovery
- [x] post-correction composite conservation validation
- [x] transactional rejection when no physical recipient exists
- [x] focused locality, conservation, and rollback coverage
- [x] GNU Fortran Debug and Release MPI qualification
- [x] all 215 serial regressions retained in MPI and non-MPI builds

This milestone removes parent-wide residual spreading from every qualified 2D
EB AMR closure path. The correction remains a conservative interface-local
closure rather than a claim of exact AMReX `MLStateRedistribute` transfer
parity.

## Embedded-wall thermal and viscous transport (`0.180.0`)

- [x] validated embedded-wall record in the shared 2D boundary set
- [x] backward-compatible adiabatic slip and impermeable defaults
- [x] physical fluid-centroid to wall-centroid normal distance
- [x] isothermal Fourier heat transfer on cut faces
- [x] no-slip Newtonian normal and tangential momentum transfer
- [x] moving-wall viscous work in total energy
- [x] exact zero mass and species wall transport
- [x] wall-length and cut-fluid-volume conservative source scaling
- [x] single-level, fixed-depth, multipatch, arbitrary-depth path reuse
- [x] serial and sparse-MPI boundary-control consensus
- [x] focused heat, traction, work, slip, locality, and rejection gates
- [x] GNU Fortran Debug and Release MPI qualification
- [x] all 215 serial regressions retained in MPI and non-MPI builds

This milestone adds the embedded diffusive flux to the established transport
right-hand side without changing its AMR or MPI stepping interfaces. Embedded
species conversion/catalytic kinetics, higher-order wall-normal gradients,
public namelist controls, and checkpoint compatibility for nondefault wall
values remain separate work.

## Public single-level embedded-wall configuration (`0.181.0`)

- [x] `embedded_wall_kind` slip/no-slip namelist control
- [x] `embedded_wall_thermal` adiabatic/isothermal namelist control
- [x] positive `embedded_wall_temperature` validation
- [x] finite three-component `embedded_wall_velocity` validation
- [x] thermal-conduction requirement for isothermal walls
- [x] viscosity requirement for no-slip walls
- [x] no-slip requirement for moving walls
- [x] transactional boundary-set configuration
- [x] public single-level EB application integration
- [x] public cut-cell heating and tangential wall-momentum regression
- [x] direct-API invalid-control rejection
- [x] active nondefault AMR application preflight rejection
- [x] GNU Fortran Debug and Release MPI qualification
- [x] all 215 serial regressions retained in MPI and non-MPI builds

Checkpoint-capable AMR application exposure remains deferred until every
nondefault embedded-wall value participates in restart compatibility. The
underlying AMR and sparse-MPI transport kernels already consume an explicitly
configured boundary set.

## Restart-safe AMR embedded-wall configuration (`0.182.0`)

- [x] one configured boundary builder shared by public EB drivers
- [x] two-level, multipatch, three-level, and arbitrary-depth wall activation
- [x] sparse-MPI public wall activation
- [x] fixed-depth checkpoint schema-2 wall records
- [x] transport enable/process/CFL restart compatibility
- [x] arbitrary-depth checkpoint fingerprint schema 3
- [x] exact wall string and tolerance-aware real-value restart matching
- [x] active isothermal-wall two-level checkpoint/restart regression
- [x] public moving no-slip/isothermal two-level AMR transport regression
- [x] wall-temperature mismatch transactional rejection
- [x] arbitrary-depth wall-fingerprint mismatch regression
- [x] GNU Fortran Debug and Release MPI qualification
- [x] all 215 serial regressions retained in MPI and non-MPI builds

The public AMR applications no longer reject transport checkpointing solely
because transport is enabled. Restart remains strict: older checkpoint schemas
and any changed wall or transport controls are rejected. Catalytic species wall
transfer and higher-order wall-normal stencils remain outside this milestone.

## EB-safe limited-linear AMR prolongation (`0.183.0`)

- [x] public reactive EB limited-linear prolongation API
- [x] component-wise MC slopes from regular coarse neighbors
- [x] exact zero-mean Cartesian child offsets per regular parent
- [x] conserved-state interpolation with EOS temperature recovery
- [x] PCM fallback for cut and topology-mismatched parents
- [x] parent-local PCM retry after an inadmissible linear candidate
- [x] exact neutral outputs after nonfinite-input rejection
- [x] linear-field child-value regression
- [x] prolongation/average-down conservation regression
- [x] explicit cut-parent fallback regression

At the `0.183.0` boundary, the kernel was qualified as a low-level alternative
to the existing PCM initializer. Namelist selection, restart fingerprinting of
that selection, and public regrid integration were deferred to later lifecycle
work.

## Fixed-depth public prolongation selection (`0.184.0`)

- [x] `pcm`/`linear` `&eb_amr` namelist selection with PCM default
- [x] one validating dispatcher for both prolongation kernels
- [x] static and dynamically regridded two-level propagation
- [x] separated sibling-patch initialization and regrid propagation
- [x] three-level coarse-to-middle and middle-to-finest propagation
- [x] public hot-wall transport case using limited-linear prolongation
- [x] invalid-method transactional rejection
- [x] linear fixed-depth checkpoint/restart preflight rejection
- [x] arbitrary-depth serial and sparse-MPI PCM-only guard

Checkpoint files do not yet identify the selected method, so linear mode is
deliberately limited to fresh fixed-depth runs. Adding the method to every
fixed-depth schema and the shared patch-tree fingerprint remains the next
lifecycle milestone before linear arbitrary-depth or restart claims.

## Restart-safe arbitrary-depth prolongation selection (`0.185.0`)

- [x] fixed-depth checkpoint schemas advanced to version 3
- [x] method stored and matched for single-patch and sibling-patch restart
- [x] method stored and matched for static and dynamic three-level restart
- [x] shared serial/sparse patch-tree fingerprint schema 4
- [x] serial arbitrary-depth initialization and iterative regrid propagation
- [x] sparse-MPI owner-local linear prolongation and direct child routing
- [x] collective sparse-MPI method validation and rank consensus
- [x] public four-level serial and 1/2/4/8-rank linear regrid qualification
- [x] two-rank linear checkpoint with four-/eight-rank restart parity
- [x] fixed-depth and arbitrary-depth method-mismatch rejection

PCM remains the public default. Limited-linear slopes are still restricted to
topology-consistent regular parent/child blocks; EB-cut parents and
inadmissible candidates keep the established PCM fallback.

## Conservative limited-linear cut-parent prolongation (`0.186.0`)

- [x] active-neighbor fluid-centroid slopes for EB-cut parents
- [x] MC limiting with one-sided embedded-boundary derivatives
- [x] fine-fluid-centroid reconstruction offsets
- [x] fine-volume-fraction-weighted zero-mean correction
- [x] component bounds from the active 3-by-3 coarse neighborhood
- [x] exact EB average-down recovery of the cut-parent conserved state
- [x] nonconstant active-child cut-parent regression
- [x] retained EOS recovery and parent-local PCM retry
- [x] retained PCM behavior for covered and topology-mismatched parents
- [x] EOS-validated order-2 StateRedist with conservative order-0 retry
- [x] dynamic three-level reactive transport/regrid qualification

Configured limited-linear initialization now remains nonconstant next to a cut
embedded boundary while preserving the cut-parent volume-weighted average.
If its sharper gradients expose an inadmissible second-order StateRedist
candidate, the whole provisional state is redistributed again with the
conservative zeroth-order kernel before commit. PCM remains the public default.
Exact AMReX EB interpolation, least-squares multidimensional gradients, and
higher-order cut-cell reconstruction are not claimed.

## Multidimensional cut-parent prolongation (`0.187.0`)

- [x] connected active coarse-cell 3-by-3 stencil
- [x] fluid-volume-centroid least-squares normal equations
- [x] full-rank two-dimensional gradient recovery
- [x] minimum-norm rank-one tangent recovery
- [x] coarse-neighbor and fine-child component envelope limiting
- [x] exact fluid-centroid tangential-affine cut-parent reproduction
- [x] retained fine-volume-weighted conservation and EOS/PCM transaction
- [x] shared fixed-depth, arbitrary-depth, serial, and sparse-MPI dispatcher

Configured limited-linear initialization now retains a genuinely
multidimensional smooth gradient at an EB-cut parent. Diagonal neighbors enter
only through an active face-connected path, so the fit does not cross a
covered corner. Rank-deficient stencils retain variation only in their
resolved tangent direction. Exact AMReX interpolation, quadratic cut-parent
reconstruction, and geometry coarsening are not claimed.

## Rank-recovering cut-parent stencil (`0.188.0`)

- [x] compact connected 3-by-3 stencil retained as the first choice
- [x] rank-deficiency detection before accepting a tangent-only gradient
- [x] bounded face-connected 5-by-5 stencil growth
- [x] disconnected and covered coarse cells excluded from the grown fit
- [x] full-rank two-dimensional recovery after a local rank-one stencil
- [x] selected-stencil component envelope limiting
- [x] exact affine reproduction through a turning fluid path
- [x] retained fine-volume-weighted conservation and EOS/PCM transaction

A cut parent no longer discards a smooth gradient component merely because
its immediate fluid neighborhood is locally one-dimensional. The stencil
grows once, remains topology connected, and falls back to a minimum-norm
rank-one fit only when the larger neighborhood is also deficient. Exact AMReX
interpolation, quadratic reconstruction, and geometry coarsening are not
claimed.

## Transactional fixed three-level parent regrid (`0.189.0`)

- [x] root-temperature plan for a replacement root-to-middle patch
- [x] deepest-first restriction of the accepted finest contribution
- [x] same-resolution middle overlap retained by the two-level regrid kernel
- [x] replacement finest plan constrained to a two-cell middle margin
- [x] configured PCM or limited-linear prolongation at both rebuilt interfaces
- [x] before/after three-level composite conservation acceptance gate
- [x] atomic publication of all fields, geometries, and patch descriptors
- [x] invalid-control rollback regression

The fixed-depth library can now relocate or resize its parent refined patch
without publishing an incomplete hierarchy. Old finest information is first
restricted into the middle level; a replacement finest patch is then
prolonged from the rebuilt middle. Public scheduling and checkpoint storage of
the moving parent topology remain the next fixed-depth lifecycle boundary.

## Public dynamic three-level parent lifecycle (`0.190.0`)

- [x] opt-in `dynamic_parent_regridding` namelist control
- [x] minimum parent-size validation for a safe nested finest plan
- [x] parent-first initialization and accepted-step regrid schedule
- [x] finest-only fallback when the scheduled parent is unchanged
- [x] one committed regrid count per parent-or-finest topology event
- [x] dynamic three-level checkpoint schema 4
- [x] stored parent-policy compatibility check
- [x] actual parent and finest descriptor reconstruction on restart
- [x] uninterrupted/checkpoint-stop/restart field parity gate

The public fixed-depth application can now move both refined rectangles while
retaining the existing fixed-parent behavior as the default. Restart rebuilds
the middle geometry from the stored parent descriptor and validates the stored
finest descriptor against that actual geometry before publishing any field.

## Arbitrary-depth outflow-boundary AMR children (`0.191.0`)

- [x] one-sided root and descendant temperature-gradient tagging
- [x] domain-inclusive recursive child planning through four levels
- [x] current-fine-state outflow exterior fill on the physical side
- [x] coarse-time interpolation on every remaining coarse/fine side
- [x] physical-side flux-register and reflux omission
- [x] serial public application topology and field gate
- [x] sparse-MPI 1/2/4/8-rank topology and field parity
- [x] retained composite EB classes, positivity, and species closure

The arbitrary-depth public lifecycle can now refine an outflow physical
boundary in serial or sparse-MPI execution. Non-outflow refined boundaries,
periodic-seam children, and physical-boundary checkpoint continuation remain
outside this milestone.

## Boundary-touching patch-tree checkpoint/restart (`0.192.0`)

- [x] four-level outflow-boundary checkpoint topology
- [x] selected-root checkpoint write and terminal marker
- [x] independent serial process restart
- [x] exact boundary reach retained at every populated level
- [x] uninterrupted/restarted composite topology parity
- [x] uninterrupted/restarted numeric field parity
- [x] two-rank checkpoint with four- and eight-rank restart
- [x] ownership-weight change across sparse restart
- [x] retained incompatible-fingerprint rejection

The qualified arbitrary-depth outflow-boundary tree can now cross a process
and MPI ownership boundary without moving away from the physical side. The
checkpoint schema remains version 4 because existing child bounds and geometry
records already describe the topology completely.

## Boundary-touching arbitrary-depth thermal transport (`0.193.0`)

- [x] public fresh boundary tree with thermal conduction enabled
- [x] recursive transport execution through four populated levels
- [x] current-fine physical-side diffusive exterior state
- [x] `r^2` child transport subcycling
- [x] physical-side diffusive register and reflux omission
- [x] serial and sparse-MPI 1/2/4/8-rank field parity
- [x] transport-active checkpoint fingerprint
- [x] independent serial transport restart parity
- [x] two-rank checkpoint to four-/eight-rank transport restart

Thermal conduction is now qualified on the arbitrary-depth outflow-boundary
topology, including process and sparse ownership changes. Viscosity, species
diffusion, and barodiffusion remain disabled in this focused public case but
retain their established interior-tree library gates.

## Boundary-touching arbitrary-depth mixture transport (`0.194.0`)

- [x] Newtonian viscosity on every populated boundary-tree level
- [x] Fourier conduction retained with the full transport combination
- [x] mixture-averaged species diffusion and correction velocity
- [x] barodiffusion and species enthalpy flux
- [x] physical-side diffusive register and reflux omission
- [x] serial and sparse-MPI 1/2/4/8-rank fresh field parity
- [x] independent serial full-transport restart parity
- [x] two-rank checkpoint to four-/eight-rank full-transport restart
- [x] schema-4 fingerprint coverage for every active transport control

The public boundary topology now qualifies the complete currently implemented
mixture-transport combination across recursive AMR, process ownership, and
checkpoint boundaries. This remains the dilute mixture-averaged model; it is
not Stefan--Maxwell, Soret, Dufour, or full PelePhysics transport parity.

## Reacting boundary-touching arbitrary-depth full physics (`0.195.0`)

- [x] elementary chemistry on both reaction half-steps at every populated level
- [x] viscosity, Fourier conduction, species diffusion, and barodiffusion
- [x] recursive hydro and transport subcycling between reaction half-steps
- [x] one transactional `R-T-H-T-R` candidate per accepted root interval
- [x] exact x-upper physical-side contact through fresh and restart lifecycles
- [x] serial and sparse-MPI 1/2/4/8-rank fresh field parity
- [x] independent serial reacting restart parity
- [x] two-rank checkpoint to four-/eight-rank reacting restart
- [x] schema-4 chemistry and transport fingerprint coverage

The public boundary tree now qualifies the selected elementary chemistry,
every currently implemented transport term, and hydro together across
arbitrary depth, process ownership, and restart. This is not arbitrary
mechanism parsing or CVODE parity.

## Restart-persistent patch-tree limiter diagnostic (`0.196.0`)

- [x] schema-2 base and schema-5 fingerprinted checkpoint envelopes
- [x] finite `[0,1]` minimum transport theta metadata validation
- [x] serial checkpoint round trip and failure rollback
- [x] selected-root sparse checkpoint write
- [x] communicator-wide restart metadata broadcast
- [x] changed-rank continuation of the pre-checkpoint minimum
- [x] retained topology, fingerprint, and field transaction boundaries

Serial and sparse-MPI patch-tree applications now report the minimum limiter
over the complete logical run, not only the post-restart suffix. This is a
diagnostic schema change; numerical state and fingerprint contents are
unchanged.

## Restart-persistent conservation diagnostic (`0.197.0`)

- [x] schema-3 base and schema-6 fingerprinted checkpoint envelopes
- [x] finite `nvar`-component initial composite-integral record
- [x] serial exact round trip and unallocated failure rollback
- [x] selected-root sparse write and communicator-wide restart broadcast
- [x] rank consensus for an explicitly supplied baseline
- [x] changed-rank continuation from the pre-checkpoint run baseline
- [x] uninterrupted/restarted conservation and limiter diagnostic comparison

Serial and sparse-MPI patch-tree applications now measure their final
composite change from the beginning of the complete logical run. Checkpoint
compatibility remains strict.

## Restart-persistent AMR operator counters (`0.198.0`)

- [x] schema-4 base and schema-7 fingerprinted checkpoint envelopes
- [x] common nonnegative chemistry/transport/hydro counter capacity
- [x] capacity larger than the currently populated tree depth
- [x] serial exact round trip and unallocated failure rollback
- [x] sparse owner-local delta reduction into global counters
- [x] selected-root sparse write and communicator-wide restart broadcast
- [x] collective rejection of rank-disagreed counter vectors
- [x] uninterrupted/restarted public counter comparison

Public serial and sparse-MPI patch-tree applications now report cumulative
per-level operator work over the complete logical run. Counter storage follows
configured maximum depth rather than current topology depth, so dynamic
coarsening and later level recreation retain earlier counts.

## Restart-persistent AMR regrid history (`0.199.0`)

- [x] schema-5 base and schema-8 fingerprinted checkpoint envelopes
- [x] cumulative successful tag/regrid evaluation count
- [x] cumulative tagged-cell count across all evaluated parent patches
- [x] serial exact round trip and neutral failure rollback
- [x] sparse rank consensus and selected-root metadata broadcast
- [x] collective rejection of rank-disagreed regrid history
- [x] uninterrupted/restarted serial diagnostic comparison
- [x] two-to-four/eight-rank sparse restart diagnostic comparison

The public adaptation diagnostics now cover the complete logical run rather
than only the current process suffix. A successful evaluation that leaves the
topology unchanged is still counted; a failed transaction is not.

## Public branching patch-tree lifecycle (`0.200.0`)

- [x] separated boundary and interior temperature-tag features
- [x] at least two leaf-visible patch identities on one populated level
- [x] four populated levels with x-upper contact on every level
- [x] complete reacting `R-T-H-T-R` physics on both branches
- [x] serial and sparse-MPI 1/2/4/8-rank fresh field parity
- [x] independent serial checkpoint/restart field parity
- [x] two-rank checkpoint continuation at four and eight ranks
- [x] retained schema-8 fingerprint and cumulative diagnostics

The public application lifecycle now qualifies a branching numerical tree,
not only a single parent-child chain. The second branch remains independent of
the physical boundary branch across regridding, ownership changes, and restart.

## Completion and toolchain contract (`0.201.0`)

- [x] recursive PeleC reference manifest with immutable commit SHAs
- [x] ordered completion workstreams with explicit exit gates
- [x] synchronized CMake, runtime-banner, README, and validation versions
- [x] Debug, Release, live-Cantera, MPI Debug, and MPI Release presets
- [x] configure-time rejection of an unusable or mixed-toolchain `mpi_f08`
- [x] explicit serial and sparse patch-tree geometry callback context
- [x] module-scope callback implementations without GNU stack trampolines
- [x] final ELF `PT_GNU_STACK` regression gate

This milestone hardens provenance, configuration, and binary security for the
existing subset. Three-dimensional flow, direct mechanism ingestion, LES,
particles/spray, scalable production I/O, and accelerators remain open on the
completion roadmap.

## Three-dimensional coordinate core (`0.202.0`)

- [x] uniform x/y/z cell-center and spacing construction
- [x] reversible z-to-x conserved and primitive state rotations
- [x] ideal-gas z-normal physical, Rusanov, and PeleC-style fluxes
- [x] passive-multispecies z-normal flux and exact species-flux closure
- [x] general-EOS reacting z-normal HLLC equal-state consistency
- [x] explicit test boundary separating the coordinate core from a 3D solver
- [ ] conservative 3D time advance and CFL selection
- [ ] 3D convergence, transport, AMR, restart, and MPI rank parity

The workstream is now in progress, not complete. No `pelef3d` executable or
three-dimensional field evolution is claimed by this milestone.

## Conservative periodic 3D Euler (`0.203.0`)

- [x] periodic x/y/z finite-volume flux divergence
- [x] summed multidimensional spectral-radius CFL selection
- [x] transactional two-stage SSPRK2 update with physical-state rejection
- [x] exact x/y/z whole-step reduction for Rusanov and PeleC-style fluxes
- [x] 16/32/64-cell-per-axis entropy-wave convergence toward first order
- [x] roundoff-level mass, three-momentum, and energy conservation
- [x] strict namelist-driven `pelef3d` application and deterministic CSV check
- [x] general-EOS reacting-state thermodynamics in 3D (`0.204.0`)
- [x] cell-local chemistry and transactional Strang splitting in 3D (`0.205.0`)
- [x] periodic regular-grid molecular transport in 3D (`0.206.0`)
- [x] 3D characteristic PLM on regular and serial static-AMR hydro (`0.210.0`)
- [ ] sparse/distributed high-order 3D AMR, dynamic topology, and EB

The regular 3D Euler update is now runnable and numerically qualified within
this boundary. It does not complete the wider three-dimensional workstream.

## General-EOS multispecies 3D Euler (`0.204.0`)

- [x] runtime NASA7 reactive state across periodic x/y/z divergence
- [x] Rusanov, HLLC, and PeleC-style general-EOS directional fluxes
- [x] summed frozen-mixture spectral-radius CFL selection
- [x] transactional SSPRK2 state and temperature recovery
- [x] exact x/y/z reduction against an independent 1D SSPRK2 reference
- [x] 6/12/24-cell-per-axis seven-species entropy-wave convergence
- [x] roundoff conservation of all Euler and species components
- [x] strict elementary/full-H2O2 public application and CSV contract
- [ ] 3D chemistry splitting and molecular transport
- [x] 3D characteristic PLM on regular and serial static-AMR hydro (`0.210.0`)
- [ ] sparse/distributed high-order 3D AMR, dynamic topology, and EB

The `reactive` name denotes the state layout and thermodynamic closure here.
No reaction source is advanced in this milestone; reacting-flow qualification
begins only when the chemistry split is added and tested separately.

## Cell-local chemistry splitting in 3D (`0.205.0`)

- [x] elementary seven-species explicit chemistry on every 3D cell
- [x] full ten-species implicit chemistry on every 3D cell
- [x] exact reduction to the established independent-cell reactor result
- [x] transactional whole-field chemistry commit with state/temperature rollback
- [x] outer `R(dt/2)-H(dt)-R(dt/2)` transaction and hydro-failure rollback
- [x] homogeneous-field reduction to two one-dimensional reactor half-steps
- [x] H/O/N elemental conservation for both mechanism paths
- [x] strict periodic Gaussian-hotspot application and independent CSV gate
- [ ] molecular transport in 3D
- [x] 3D characteristic PLM on regular and serial static-AMR hydro (`0.210.0`)
- [ ] sparse/distributed high-order 3D AMR, dynamic topology, and EB

The regular 3D path now qualifies a reacting-flow source split for the two
built-in H2/O2 mechanisms. This does not establish arbitrary mechanism
ingestion, external ignition validation, molecular transport, or whole-program
PeleC parity.

## Regular-grid molecular transport in 3D (`0.206.0`)

- [x] periodic x/y/z diffusive finite-volume divergence
- [x] full 3-by-3 Newtonian velocity-gradient stress with Stokes hypothesis
- [x] Fourier conduction and species-enthalpy energy flux
- [x] mixture-averaged species diffusion and barodiffusion
- [x] correction velocity with one conservative upper face per direction
- [x] six-face outgoing species-positivity limiter
- [x] summed three-dimensional parabolic stable-step selection
- [x] transactional transport SSPRK2 state and temperature recovery
- [x] outer `R-T-H-T-R` full-physics transaction and late-failure rollback
- [x] elementary/full-H2O2 x/y/z reduction against established 2D transport
- [x] second-order diagonal viscous-wave convergence
- [x] independent thermal and species smoothing with periodic conservation
- [x] public full-H2O2 transport/control driver and CSV comparison
- [x] 3D characteristic PLM on regular and serial static-AMR hydro (`0.210.0`)
- [ ] sparse/distributed high-order 3D AMR, dynamic topology, and EB

This closes molecular transport only for the serial, uniform, single-level,
periodic 3D path. Soret/Dufour, multicomponent diffusion, arbitrary transport
data ingestion, nonperiodic walls, and production PeleC transport parity are
not established.

## Static two-level reactive 3D AMR hydro (`0.207.0`)

- [x] validated strictly interior rectangular patch descriptor
- [x] conservative PCM prolongation and volume-average restriction
- [x] composite all-component integral and average-down operations
- [x] SSPRK2 time-averaged coarse x/y/z face-flux publication
- [x] integer-ratio fine-level subcycling
- [x] linearly time-interpolated coarse boundary state on six fine faces
- [x] face-area and substep-time averaged fine interface fluxes
- [x] six-face coarse/fine reflux followed by average-down
- [x] coarse/fine hyperbolic CFL selection in root-step units
- [x] transactional state and temperature rollback on late rejection
- [x] uniform-state invariance and nonuniform composite conservation
- [x] public full-H2O2 application with independent coarse/fine CSV check
- [ ] chemistry and molecular-transport coupling across AMR levels
- [x] characteristic PLM and limited coarse ghost prolongation (`0.210.0`)
- [ ] CTU/PPM reconstruction and dynamic regridding
- [x] versioned transactional two-level checkpoint/restart (`0.208.0`)
- [x] replicated-state MPI x-slab decomposition and rank parity (`0.209.0`)
- [ ] sparse MPI hierarchy storage and 3D EB

This milestone begins the 3D AMR workstream with one serial, static,
strictly interior patch over a periodic coarse grid. It isolates hydrodynamic
flux synchronization before introducing sources, diffusion, topology change,
or distributed ownership.

## Static 3D AMR checkpoint/restart (`0.208.0`)

- [x] self-identifying checkpoint magic and explicit schema version
- [x] complete NASA7 species-order and coefficient fingerprint
- [x] solver, reconstruction, limiter, boundary, mesh, patch, ratio, CFL, and physics fingerprint
- [x] both conserved levels and recovered-temperature payloads
- [x] physical time, cumulative step, conservation baseline, and reflux history
- [x] write-time physicality, synchronization, and metadata rejection
- [x] private candidate load with commit only after complete validation
- [x] truncated, trailing, incompatible, and output-path collision rejection
- [x] configured interval checkpoint and intentional stop-after-checkpoint
- [x] public full-H2O2 continuous/restart byte-exact coarse/fine CSV parity
- [ ] dynamic or multipatch 3D topology checkpoint
- [x] selected-root two-to-four/eight-rank restart (`0.209.0`)
- [ ] distributed checkpoint layout

This milestone closes restart only for the qualified serial, static,
single-patch, hydro-only 3D AMR hierarchy. It does not imply MPI ownership,
chemistry/transport AMR state, EB metadata, portable cross-schema conversion,
or scalable production I/O.

## Distributed-slab static 3D AMR hydro (`0.209.0`)

- [x] deterministic contiguous coarse/fine x-slab ownership
- [x] one or more coarse and fine x planes per participating rank
- [x] owner-only multidimensional CFL scans and communicator maximum
- [x] owner-only x/y/z face evaluation and SSPRK2 cell updates
- [x] collectively reconstructed face and state candidates
- [x] serial-identical six-face reflux and conservative average-down
- [x] transactional rejection of invalid solver, shape, root, and metadata
- [x] bitwise NASA7 name/range/coefficient consensus across ranks
- [x] selected-root initialization, checkpoint read/write, and CSV output
- [x] serial and 1/2/4/8-rank byte-exact public coarse/fine parity
- [x] two-rank checkpoint continuation at four and eight ranks
- [ ] rank-local persistent hierarchy storage and halo exchange
- [ ] distributed/scalable checkpoint and field output
- [ ] AMR chemistry, molecular transport, distributed PLM, regridding, and EB

This milestone closes the completion roadmap's first distributed arithmetic
gate for the qualified static 3D AMR hydro case. It deliberately keeps both
levels replicated, so it is reproducibility evidence rather than a production
scaling claim.

## Characteristic PLM on regular and static-AMR 3D hydro (`0.210.0`)

- [x] frozen-composition characteristic slopes in x, y, and z
- [x] MC and minmod limiter configuration and rejection gates
- [x] common density, pressure, and species positivity scaling
- [x] transactional PLM SSPRK2 state, temperature, and mean-flux publication
- [x] exact y/z rotational symmetry against the x-normal PLM update
- [x] 6/12/24-cell diagonal entropy-wave orders `2.02` and `2.46`
- [x] public full-H2O2 regular-grid PLM CSV contract
- [x] characteristic PLM on both serial static AMR levels
- [x] two-layer, time-interpolated, limited-linear coarse ghost fill
- [x] PLM flux-register reflux and roundoff composite conservation
- [x] fine density L1 reduction from `1.2481e-3` (PCM) to `4.3313e-4`
- [x] schema-2 reconstruction/limiter fingerprint and exact PLM restart
- [ ] distributed AMR characteristic PLM and sparse halo exchange
- [ ] CTU transverse prediction, PPM, physical boundaries, and 3D EB

This milestone closes the selected high-order method-of-lines PLM boundary for
the periodic regular grid and serial static hierarchy. It does not upgrade the
PCM-only replicated MPI kernel or claim PeleC's multidimensional Godunov/PPM
algorithm.

## Static axis-plane reactive 3D EB hydro (`0.211.0`)

- [x] exact x/y/z planar cell volume fractions and classifications
- [x] normalized cell and all tangential open-face centroids
- [x] EB area, physical centroid, solid-to-fluid normal, and normal integral
- [x] validator-enforced zero-based face bounds and metric finiteness/ranges
- [x] cellwise Cartesian-aperture/EB-normal discrete geometric identity
- [x] explicit rejection of an interior grid-face-aligned plane
- [x] general-EOS impermeable slip-wall pressure traction
- [x] aperture-weighted x/y/z conservative flux divergence
- [x] PCM Rusanov, HLLC, and PeleC-style open-face flux construction
- [x] zero-gradient physical-domain exterior state
- [x] transactional active-cell Euler update and temperature recovery
- [x] exact x/y/z uniform tangential-flow invariance
- [x] roundoff nonuniform fluid-volume conservation and species closure
- [ ] oblique or curved 3D geometry
- [ ] grid-face-aligned interior EB ownership
- [ ] small-cell state/flux redistribution and EB-aware stable timestep
- [ ] characteristic reconstruction, CTU/PPM, chemistry, and transport
- [ ] AMR/reflux/regridding, checkpoint/restart, and MPI ownership

This milestone fixes the first 3D EB data, sign, indexing, wall-pressure, and
transaction conventions on one serial Cartesian level. It is not yet a
production cut-cell solver: without redistribution the stable interval is
limited by the smallest positive volume fraction, and no public 3D EB driver
or field-output contract is claimed.

## Conservative 3D EB FluxRedist (`0.212.0`)

- [x] validated four-dimensional residual and output extent contract
- [x] exact regular-cell identity and covered-cell zero residual
- [x] six open-face active-neighbor discovery without boundary overreach
- [x] volume-weighted cut-cell neighborhood residual
- [x] conservative return of removed excess to active neighbors
- [x] exact uniform active-residual preservation
- [x] x/y/z one-neighbor analytical cut and receiver values at `kappa=0.05`
- [x] both transverse directions and overlapping-neighborhood conservation
- [x] standalone transactional NASA7 state/temperature update
- [x] end-to-end PCM wall/divergence/FluxRedist/EOS composition
- [x] raw negative-density versus physical redistributed full-grid CFL gate
- [x] exact x/y/z uniform tangent-flow invariance after redistribution
- [x] roundoff nonuniform component conservation and species closure
- [x] invalid geometry, nonfinite residual, solver, timestep, and state gates
- [ ] weighted zeroth-/higher-order 3D StateRedist parity
- [ ] public EB-aware CFL selection and multi-step application
- [ ] general geometry and high-order EB reconstruction
- [ ] chemistry/transport, AMR/reflux/regridding, restart, and MPI ownership

This milestone removes the direct `1/kappa` failure in the qualified planar
regression without claiming that all cut-cell stability or accuracy concerns
are closed. FluxRedist remains first order, and the raw Euler route remains a
diagnostic reference rather than the preferred small-cell path.

## Zeroth-order planar 3D EB StateRedist (`0.213.0`)

- [x] distinct provisional-state redistribution API beside FluxRedist
- [x] configurable target volume fraction in `(0,1]`, defaulting to `0.5`
- [x] aperture-difference-normal selection of one regular axis receiver
- [x] explicit rejection of non-axis normals and missing regular receivers
- [x] `nrs` overlap counts and complementary self/neighbor partitions
- [x] volume-weighted zeroth-order neighborhood gather and scatter
- [x] exact uniform active-state preservation and covered standalone zero
- [x] fluid-volume conservation for every arbitrary conserved component
- [x] x/y/z closed-form cut and receiver values at `kappa=0.05`
- [x] transactional caller-residual StateRedist/NASA7 recovery
- [x] common raw/FluxRedist/StateRedist hydro residual construction
- [x] end-to-end PCM wall/divergence/StateRedist/EOS composition
- [x] physical StateRedist result at a raw-negative full-grid CFL step
- [x] nonuniform component conservation, species closure, and covered identity
- [x] invalid target, state, solver, timestep, and EOS rollback gates
- [ ] second- or higher-order 3D StateRedist
- [ ] oblique, curved, or general multi-neighbor 3D EB geometry
- [ ] public EB-aware CFL selection and multi-step application
- [ ] chemistry/transport, AMR/reflux/regridding, restart, and MPI ownership

The qualified path is deliberately order zero and exact-axis-plane only. It
closes the first weighted-state stabilization contract without turning the 3D
EB library regression into a public solver or claiming field parity with the
full AMReX StateRedist implementation.

## Public stabilized planar 3D EB application (`0.214.0`)

- [x] active-cell general-EOS full-grid CFL selector
- [x] analytical three-direction spectral-rate unit gate
- [x] covered-cell exclusion and all-covered/invalid-active-state rejection
- [x] dedicated validated `&reactive_eb_3d` input record
- [x] exact plane, regular-receiver, active-target, solver, and finite gates
- [x] public rejection of the raw unstabilized route
- [x] installed `pelef_reactive_eb_3d` executable
- [x] repeated transactional forward-Euler steps with final-time clipping
- [x] selectable public FluxRedist or order-zero StateRedist
- [x] fluid-volume diagnostics and active general-EOS extrema
- [x] wall-normal momentum separated from conserved fluid invariants
- [x] x-fastest EB geometry, primitive, and species CSV contract
- [x] independent 480-row geometry/EOS/composition/invariant checker
- [x] two-step public StateRedist and FluxRedist density-sheet regressions
- [ ] characteristic EB reconstruction or higher-order time integration
- [ ] chemistry, molecular transport, or configurable mechanism/composition
- [ ] general geometry, AMR/reflux/regridding, restart, and MPI ownership

The public application closes the single-level planar hydro workflow, not the
3D EB workstream. It is a deterministic serial validation path with explicit
method boundaries; production coupled physics and distributed hierarchy work
remain open.

## Optional chemistry in the planar 3D EB application (`0.215.0`)

- [x] optional active-cell mask on the 3D cell-local chemistry wrapper
- [x] transactional per-plane chemistry candidate and temperature recovery
- [x] transactional `R(dt/2)-H(dt)-R(dt/2)` EB composition
- [x] exact covered-cell and masked-off active-cell state/temperature identity
- [x] inert StateRedist hydro parity through the split wrapper
- [x] invalid mask, solver, and chemistry-tolerance rollback gates
- [x] active species change with mass, energy, and tangential-momentum gates
- [x] H/O/N elemental totals and species-closure/physicality gates
- [x] validated chemistry namelist controls and frozen elementary mechanism
- [x] separate public inert/reacting StateRedist cases
- [x] byte-exact inert output against the frozen 0.214 StateRedist baseline
- [x] independent 480-row geometry, EOS, covered-row, and chemistry checker
- [ ] configurable or production-stiff mechanisms and higher-order EB chemistry
- [ ] molecular transport, general geometry, AMR/reflux/regridding, restart,
  and MPI chemistry

This increment qualifies optional elementary chemistry only on the existing
serial, single-level, inviscid, exact-axis-plane PCM EB boundary. Chemistry is
active on fluid cells selected by the mask, while the embedded wall's normal
momentum remains a physical exchange rather than a conserved fluid component.
The public cases use order-zero StateRedist and the fixed seven-species
H2/O2/N2 mechanism; broader coupled-physics and distributed boundaries remain
open.

## Molecular transport in the planar 3D EB application (`0.216.0`)

- [x] active-only primitive recovery and covered-storage-independent timestep
- [x] Newtonian viscosity and Fourier conduction on open interior faces
- [x] mixture-averaged species diffusion, barodiffusion, correction velocity,
  and species-enthalpy flux
- [x] aperture-weighted divergence divided by physical fluid volume
- [x] exact zero physical-domain and adiabatic slip/impermeable EB transport
  flux
- [x] physical species-inventory outflow limiter with energy-flux consistency
- [x] transactional SSPRK2 transport with StateRedist at both Euler stages
- [x] `R-T-H-T-R` whole-step transaction and disabled-path arithmetic parity
- [x] minimum hydro/transport public step and validated process controls
- [x] exact covered state/temperature identity and invalid-input rollback
- [x] pre-assignment output-shape rejection and uniform barodiffusion dependency
- [x] physical-inventory-scaled public conservation diagnostics
- [x] x/y/z unit conservation, limiter, boundary, and split gates
- [x] public control, transport-only, and transport-plus-chemistry regressions
- [x] independent 480-row geometry, EOS, conservation, response, schedule, log,
  covered-row, and full-field hash checker
- [ ] isothermal/no-slip/catalytic or prescribed-flux 3D embedded walls
- [ ] FluxRedist transport, general geometry, AMR/reflux/regridding, restart,
  and MPI transport

The frozen public cases advance one clipped `2.0e-9 s` step, before their
zero-gradient hydro boundary becomes part of the invariant comparison.
Transport changes the active solution by `1.0489647406832604e-1`; chemistry
adds a `4.043407673691575e-8` species response. The effective relative
conservation and elemental errors are `6.867e-15` and `1.112e-14`, with maximum
diffusivity `9.7180419287635654e-4` and limiter factor `1.0`.

This qualifies the fixed elementary mixture on one serial, single-level,
exact-axis-plane StateRedist mesh. It does not claim long-time physical-
boundary validation, general reaction-aware stability, higher-order EB
transport accuracy, configurable mechanisms, or external PeleC field parity.

## Transactional planar 3D EB checkpoint/restart (`0.217.0`)

- [x] self-identifying formatted checkpoint magic, schema, and end marker
- [x] complete NASA7 thermo, elementary-reaction, and gas-transport records
- [x] immutable mesh, planar geometry, physics, and numerical fingerprints
- [x] full conserved state, recovered temperature, and x-fastest storage
- [x] initial conserved, L1, and H/O/N inventory continuation
- [x] cumulative minimum-step, maximum-diffusivity, and minimum-theta history
- [x] private candidate parsing and publish-only-after-complete validation
- [x] incompatible, truncated, and temperature-inconsistent rollback gates
- [x] default-off checkpoint cadence, stop-after-write, and restart controls
- [x] lexical output/checkpoint and output/restart alias rejection
- [x] checkpoint publication only after a committed whole split step
- [x] coupled chemistry-plus-transport stop/restart application regression
- [x] byte-exact uninterrupted/restarted final CSV and diagnostic parity
- [ ] MPI or distributed checkpoint layout
- [ ] AMR topology, reflux, or regridding state
- [ ] general 3D EB geometry or checkpoint schema conversion

This increment closes process-boundary continuation for the existing serial,
single-level, exact-axis-plane StateRedist workflow. Restart permits changed
output/checkpoint paths, cadence, and extended final step/time limits, but it
rejects changes to the model that advances the stored state. The formatted
single-file layout is deterministic validation I/O, not a scalable production
checkpoint format.
The portable path gate normalizes `/`, `.`, and reducible `..` components; it
does not claim symlink, hardlink, mount, or current-directory identity checks.

## Pinned Cantera YAML ingestion (`0.218.0`)

- [x] repository-pinned Cantera `h2o2.yaml` with exact SHA-256 provenance
- [x] mandatory explicit `ohmech` ideal-gas phase selection
- [x] ten ordered species with element composition and molecular weight
- [x] NASA7 low/high intervals at a validated 101325 Pa reference pressure
- [x] gas Lennard-Jones transport in the Fortran database units
- [x] elementary, named-collider three-body, and Troe falloff conversion
- [x] all 29 one-based source indices and six duplicate flags retained
- [x] element-balance, finite-value, identifier, ordering, and size validation
- [x] deterministic normalized JSON and generated Fortran cleanliness gates
- [x] generated thermo and transport replace the former hand-maintained tables
- [x] live Cantera trajectory plus regular 1D/2D full-H2/O2 parity retained
- [ ] runtime mechanism loading or application-level arbitrary dispatch
- [ ] CHEMKIN, NASA9, PLOG, Chebyshev, SRI, Lindemann, or custom-order support
- [ ] detailed hydrocarbon mechanisms or production CVODE-class integration

The qualified route is intentionally a build-time source pipeline. The
committed bundle records both the YAML metadata version and the Cantera API
version used for ingestion, together with the source hash, phase, and source
and target unit systems. It is not a claim that an arbitrary Cantera mechanism
can be selected at runtime.

## Selected normalized mechanism build (`0.219.0`)

- [x] optional `PELEF_MECHANISM_BUNDLE` configure-time selection
- [x] JSON-declared module and loader names read by CMake without fixed H2/O2
  identifiers
- [x] deterministic build-tree Fortran generation and separate static library
- [x] optional source-YAML basename and SHA-256 gate at configure and build
  time
- [x] automatic CMake reconfiguration after selected bundle or source changes
- [x] configured probe for thermo, reaction, primitive transport, rate, and
  Jacobian interfaces
- [x] finite-value, ordering, geometry, and common-thermo-interval checks
- [x] compiled two-species/one-reaction fixture in every test configuration
- [x] positive and negative source-hash contract tests
- [x] clean tests-disabled full-H2/O2 configure, build, run, and install gate
- [x] energy-constrained reduced-Jacobian directional finite-difference gate
- [x] compiled 128-character equation and Fortran namespace rejection gates
- [x] configure-time selected-bundle dispatch in an isolated regular serial 1D
  reacting-flow application (`0.223.0`)
- [x] configure-time selected-bundle dispatch in an isolated periodic regular
  serial 3D reacting-flow application (`0.224.0`)
- [x] configure-time selected-bundle dispatch in an isolated regular serial
  2D reacting-flow application with physical boundaries (`0.225.0`)
- [x] configure-time selected-bundle dispatch in the serial reactive 1D AMR
  application (`0.226.0`)
- [x] configure-time selected-bundle dispatch in the serial single-level
  reactive 2D EB application (`0.227.0`)
- [x] configure-time selected-bundle dispatch in the bounded regular MPI 1D
  verification application (`0.228.0`)
- [ ] selected-bundle dispatch in EB AMR/3D or MPI AMR/EB applications
- [ ] runtime YAML loading, mechanisms above 32 species, or new rate families

The optional probe is a build and ABI qualification boundary. It does not
claim that a successfully compiled mechanism is physically appropriate, nor
does it replace the applications' existing elementary/full-H2O2 selection.

## Selected normalized mechanism 0D reactor (`0.220.0`--`0.222.0`)

- [x] configure-time generated and installed `pelef0d_selected` executable
- [x] isolated reusable runtime modules without fixed generated-mechanism
  symbols
- [x] exact species-name to bundle-order mole-fraction mapping
- [x] finite scalar, strict composition-tail, duplicate, unknown-species, and
  common-thermo-interval rejection
- [x] lexical `.`, repeated-separator, and reducible-`..` input/output alias
  rejection
- [x] generator-side elemental and molecular-mass reaction-balance rejection
- [x] runtime reaction structure and molecular-mass balance checks before CSV
  creation
- [x] adaptive implicit constant-volume integration for 2--32 species
- [x] bounded initial step/output cadence with exact short final-fragment
  scheduling
- [x] schema-1 complete thermo/transport CMake contract and explicit rejection
  of generator-only legacy kinetics bundles
- [x] dynamic bundle-order mass-fraction and molar-production-rate CSV columns
  (`wdot_*` in `kmol/(m^3 s)`)
- [x] independent two-species conservation and activity checker
- [x] byte-identical repeated fixture trajectory
- [x] byte-identical selected versus fixed full-H2/O2 trajectory
- [x] selected executable installation and non-executable-stack audit
- [ ] runtime YAML/JSON loading or mechanism switching in one executable
- [x] selected-mechanism regular serial 1D chemistry and molecular-transport
  dispatch (`0.223.0`)
- [x] selected-mechanism periodic regular serial 3D chemistry and
  molecular-transport dispatch (`0.224.0`)
- [x] selected-mechanism regular serial 2D chemistry, molecular transport,
  CTU, and physical-boundary dispatch (`0.225.0`)
- [x] selected-mechanism serial reactive 1D AMR dispatch across two-level,
  multilevel, and multipatch modes (`0.226.0`)
- [x] selected-mechanism serial single-level reactive 2D EB dispatch
  (`0.227.0`)
- [x] selected-mechanism bounded regular MPI 1D dispatch (`0.228.0`)
- [ ] selected-mechanism EB AMR/3D or MPI AMR/EB dispatch
- [x] optional selected-0D SUNDIALS 7.2.0 CVODE BDF integration with pinned
  full-H2/O2 structural and Cantera validation
- [x] optional exact-version, double-precision, static SUNDIALS
  official-Fortran dependency gate, isolated from the native selected runtime
  and fixed applications, with non-double-precision rejection tests
- [x] CVODE BDF serial-vector and dense-matrix/linear-solver lifecycle
- [x] arbitrary dependent-species energy-constrained reduced Jacobian with
  centered directional finite-difference coverage
- [x] private reduced state with largest-initial-species closure and
  recoverable invalid-trial rejection
- [x] caller-state rollback on integration failure, failed-handle quarantine,
  assigned-handle reinitialization rejection, and idempotent finalization
- [x] separate native and cumulative CVODE internal-step limits, including
  pre-call exhaustion rejection, plus exact short-final scheduler behavior
- [x] deterministic two-species CVODE fixture and full 10-species/29-reaction
  H2/O2 conservation trajectory
- [x] pinned full-H2/O2 temperature, pressure, species, and instantaneous-rate
  comparison with Cantera 3.2.0
- [x] tests-disabled selected build/install and static dependency audit
- [x] opaque generation-checked handle API with 64 live serial contexts,
  per-context SUNDIALS resources, budgets, statistics, and failure state
- [x] C-interoperable callback tokens resolved through
  `FCVodeSetUserData`; no Fortran mechanism object crosses the C ABI
- [x] exact standalone/interleaved trajectory and solver-statistics identity
  for two states with different density, temperature, and closure species
- [x] peer failure/finalization isolation, stale-copy rejection after slot
  reuse, reverse-order cleanup, double finalization, and capacity exhaustion
- [x] dependency-ordered terminal cleanup with the first failing vendor
  destroy routine reported before the public handle is invalidated
- [ ] concurrent, reentrant, or thread-safe calls, sparse linear solvers,
  production performance, or detailed-fuel validation

The selected reactor consumes the generated loaders and defaults to the shared
generic elementary-reaction backward-Euler path. An optional dependency-isolated
runtime dispatches the same input and CSV contract to CVODE BDF. Its callbacks
use the generic RHS and energy-constrained semi-analytic Jacobian; the generated
production kernel remains limited to reported rates. The current SUNDIALS
adapter holds up to 64 process-local contexts in a private registry. Public
handles carry only a slot and generation; callbacks receive a stable
C-interoperable token and validate both values before resolving the Fortran
state. Context calls may be interleaved sequentially, but the registry is not
thread safe or reentrant. Successful structural and parity gates do not
establish physical suitability of a newly supplied mechanism.

## Selected normalized mechanism regular 1D reacting flow (`0.223.0`)

- [x] configured and installed `pelef_reactive_1d_selected` executable
- [x] separate selected-reactive runtime that cannot link or collide with the
  committed generated full-H2/O2 module
- [x] independent `gas_transport_mod` type/validation layer shared by fixed
  and selected primitive transport records
- [x] exact, nonduplicate species-name composition mapping into bundle order
- [x] finite active/unused composition-tail validation and safe normalization
- [x] generated NASA7, reaction, and transport loading plus startup structural,
  molecular-mass-balance, and common-temperature-range checks
- [x] existing regular 1D conservative hydro, native chemistry, viscosity,
  conduction, mixture diffusion, barodiffusion, correction velocity, and
  species-enthalpy transport
- [x] independent two-species activity, closure, conservation, uniformity, and
  deterministic-output checks
- [x] byte-identical selected/fixed ten-species, 29-reaction full-H2/O2 CSV
- [x] unknown-species and fixed-application `selected` startup rejection with
  no output artifact
- [x] tests-disabled selected build, install, stack, dependency, and executable
  identity audit
- [ ] runtime mechanism loading, selected CVODE use inside the flow solver, or
  mechanisms above 32 species
- [x] later regular 2D/3D, serial AMR 1D, serial EB 2D, and bounded regular
  MPI 1D selected dispatch (`0.224.0`--`0.228.0`)
- [ ] selected EB AMR/3D, MPI AMR/EB, thread-safe/performance, detailed-fuel,
  or external PeleC/experiment qualification

This milestone closes one application-level configure-time dispatch boundary.
The ordinary fixed executable still rejects `chemistry_model = "selected"`;
the selected executable is generated only when a normalized schema-1 bundle is
configured. It deliberately rejects AMR and checkpoint/restart controls rather
than silently ignoring them. Its result proves deterministic reuse of the
existing serial regular-grid 1D algorithm, not general PeleC completion.

## Selected normalized mechanism regular 3D reacting flow (`0.224.0`)

- [x] configured and installed `pelef_reactive_3d_selected` executable
- [x] shared exact-name composition validation and safe normalization across
  selected 1D and 3D configuration front ends
- [x] selected mechanism thermo/reaction/mass/transport/common-range startup
  validation consolidated in a reusable dependency-isolated module
- [x] separate selected 3D module directory and link graph without committed
  generated elementary or full-H2/O2 symbols
- [x] shared fixed/selected regular 3D application driver receiving validated
  species, reactions, transport, and bundle-order base composition
- [x] periodic PCM or characteristic-PLM SSPRK2 hydro, native chemistry, and
  complete existing regular-3D molecular-transport composition
- [x] independently ordered two-species 4-by-4-by-4 chemistry fixture with
  closure, physicality, uniformity, activity, and repeat-byte gates
- [x] full-H2/O2 selected transport hotspot passes the independent fixed
  control/transport response checker
- [x] byte-identical selected/fixed full-H2/O2 transport-hotspot CSV
- [x] active full-H2/O2 chemistry through a nonuniform characteristic-PLM
  hotspot with byte-identical selected/fixed CSV
- [x] unknown-species and fixed-application `selected` startup rejection with
  no output artifact
- [x] later regular 2D physical boundaries, serial AMR 1D, serial EB 2D, and
  bounded regular MPI 1D selected dispatch (`0.225.0`--`0.228.0`)
- [ ] 3D physical boundaries, selected EB AMR/3D, MPI AMR/EB, runtime loading,
  CFD CVODE, mechanisms above 32 species, thread safety/performance, or
  detailed fuels

This milestone closes the periodic single-level regular-3D configure-time
dispatch boundary. One full-H2/O2 parity case resolves nonzero conduction,
species/barodiffusion, and flow response; a second executes active chemistry
through characteristic PLM, while the smaller independent bundle proves a
different generated namespace and ordering. These are deterministic
implementation gates, not external PeleC field parity or physical mechanism
validation.

## Selected normalized mechanism regular 2D reacting flow (`0.225.0`)

- [x] configured and installed `pelef_reactive_2d_selected` executable
- [x] shared fixed/selected application driver for simulation, deterministic
  CSV publication, invariants, extrema, and diagnostics
- [x] private selected 2D runtime without either committed generated mechanism
- [x] exact-name composition resolution for 2--32 species
- [x] 32-slot prescribed wall-flux storage with actual-species zero-tail
  validation after the selected bundle is loaded
- [x] regular PCM, characteristic PLM/PPM, optional CTU, native chemistry,
  molecular transport, and configured periodic/physical boundary paths
- [x] independently ordered two-species active-chemistry fixture with positive
  state, closure, uniformity, and repeat-byte gates
- [x] byte-identical fixed/selected full-H2/O2 uniform chemistry output
- [x] byte-identical fixed/selected full-H2/O2 characteristic-PLM, CTU, active
  chemistry output
- [x] byte-identical fixed/selected full-H2/O2 prescribed-species-wall
  transport output
- [x] unknown-species, duplicate, nonfinite-tail, composition-wave endpoint,
  wall-flux-tail, and fixed-application startup rejection
- [x] selected serial single-level reactive 2D EB dispatch (`0.227.0`)
- [x] later bounded regular MPI 1D selected dispatch (`0.228.0`)
- [ ] selected EB AMR/3D, MPI AMR/EB, runtime loading, CFD CVODE, mechanisms
  above 32 species, thread safety/performance, or detailed-fuel qualification

The selected 2D target links the selected 1D runtime and privately compiles
only the regular multidimensional configuration, boundary, transport,
evolution, and application modules. The selected 3D runtime now builds on this
2D target instead of compiling a second copy of the same modules. Complete
schema-1 bundles explicitly select the generic adaptive `explicit` or
`implicit` path; fixture and full-H2/O2 bundles pin those policies,
respectively. Neither path calls a mechanism-specific reactor or infers the
selected policy from species count.

## Selected normalized mechanism serial AMR 1D reacting flow (`0.226.0`)

- [x] configured and installed `pelef_amr_reactive_1d_selected` executable
- [x] shared fixed/selected application driver for AMR mode dispatch,
  simulation, deterministic composite CSV, conservation, and diagnostics
- [x] private selected AMR runtime without either committed mechanism
- [x] optional bundle-order composition threaded through two-level,
  arbitrary-depth, and dynamic multipatch root initialization
- [x] existing PCM, PLM, PPM, characteristic PPM, native chemistry, molecular
  transport, subcycling, reflux, average-down, and regridding algorithms
- [x] independently ordered two-species active chemistry/transport fixture
  with refinement, closure, composite coverage, and repeat-byte gates
- [x] byte-identical fixed/selected full-H2/O2 two-level active
  chemistry/transport output
- [x] byte-identical fixed/selected full-H2/O2 three-level characteristic-PPM
  output
- [x] byte-identical fixed/selected full-H2/O2 three-patch molecular-transport
  output
- [x] three-level boundary-touching outflow with characteristic PPM and hybrid
  WENO7-Z, byte-identical between fixed and selected front ends
- [x] explicit bundle chemistry-integrator policy threaded through selected
  regular 1D/2D/3D and AMR 1D without species-count inference
- [x] unknown-species, non-AMR, checkpoint/restart, input/output alias,
  entropy-temperature-range, and fixed-application startup rejection without
  an output artifact
- [x] later bounded regular MPI 1D selected dispatch (`0.228.0`)
- [ ] selected AMR checkpoint/restart, EB AMR/3D, MPI AMR/EB, runtime loading,
  CFD CVODE, mechanisms above 32 species, thread safety/performance, or
  detailed fuels

The target links the selected regular-1D runtime and privately compiles only
the generic hierarchy, regrid, reacting AMR, and application modules. Fixed
callers omit the new optional composition argument and preserve their prior
byte-exact public results. This dispatch evidence does not qualify a selected
mechanism physically or establish external PeleC field parity.

## Selected normalized mechanism serial reactive EB 2D flow (`0.227.0`)

- [x] configured and installed `pelef_reactive_eb_2d_selected` executable
- [x] shared fixed/selected application driver for geometry and boundary
  construction, simulation, deterministic EB CSV, volume-weighted invariants,
  and diagnostics
- [x] private selected EB runtime without either committed mechanism
- [x] optional bundle-order composition used by regular initialization and
  configured physical-boundary states
- [x] generated explicit/implicit chemistry policy used by both active-cell
  Strang half-steps with complete rollback on policy rejection
- [x] existing plane/circle geometry, PCM/characteristic-PLM hydro,
  FluxRedist/StateRedist controls, native chemistry, molecular transport, and
  isothermal/no-slip embedded-wall paths
- [x] explicit outflow-only outer-domain boundary constraint retained
- [x] independently ordered two-species active chemistry/transport plane
  fixture with 8 regular, 4 cut, and 4 covered cells and repeat-byte output
- [x] byte-identical fixed/selected full-H2/O2 plane chemistry output
- [x] byte-identical fixed/selected full-H2/O2 characteristic-PLM circle
  transport output with an isothermal moving no-slip embedded wall
- [x] unknown-species, embedded-wall-temperature, input/output alias, and
  fixed-application selected-model rejection before output creation
- [x] later bounded regular MPI 1D selected dispatch (`0.228.0`)
- [ ] selected EB AMR/3D, checkpoint/restart, MPI AMR/EB, runtime loading,
  CFD CVODE, general outer physical boundaries, mechanisms above 32 species,
  thread safety/performance, or detailed fuels

The selected target links the selected regular-2D runtime and privately
compiles the generic EB geometry, reconstruction, hydro, transport, driver,
and application modules. The evidence closes configure-time dispatch only for
the established serial single-level 2D EB boundary. It does not qualify a new
mechanism physically or establish external PeleC field parity.

## Selected normalized mechanism regular MPI 1D flow (`0.228.0`)

- [x] configured and installed `pelef_mpi_reactive_1d_selected` executable
- [x] shared fixed/selected driver for the uneven 19-cell decomposition,
  deterministic initialization, adaptive `R-T-H-T-R` loop, collective
  invariants, ordered gather, dynamic-species CSV, and diagnostics
- [x] private selected MPI runtime without `pelef_core`, fixed database
  loaders, or either committed generated mechanism
- [x] generated explicit/implicit chemistry policy threaded through both
  distributed chemistry half-steps and adaptive retries
- [x] communicator-wide policy consensus plus collective invalid or
  rank-disagreement failure with exact caller state and temperature rollback
  through a nonzero adaptive Strang call on 1/2/4 ranks
- [x] independently ordered two-species explicit fixture with active chemistry,
  positive state, closure, conservation, and byte-identical 1/2/4-rank output
- [x] selected implicit full-H2/O2 output byte-identical across 1/2/4 ranks and
  to the fixed MPI executable
- [x] refactored fixed 1/2/4-rank output byte-identical to the frozen
  pre-refactor one-rank CSV
- [x] dynamic-species MPI comparator and generator namespace/template gates
- [x] sanitized tests-disabled Release build/install with 36 build/install
  byte-identical ELF executables, non-executable stacks, no RPATH/RUNPATH, and
  one coherent OpenMPI dependency set
- [x] configure-time wrapper/launcher co-location and post-link unresolved or
  mixed-MPI-ABI rejection
- [x] exact ordered H2/H and pinned ten-species full-H2/O2 initialization
  profiles; other selected MPI bundle orders are rejected
- [ ] general input-driven selected MPI flow, physical boundaries, selected
  MPI AMR/EB, checkpoint/restart, runtime loading, CFD CVODE, mechanisms above
  32 species, thread safety/scaling, detailed fuels, or external validation

The selected target links the selected regular-1D runtime and privately
compiles the generic MPI domain, reactive transport, reactive advance, and
shared application modules. Its zero-or-one-output-filename interface and
bundle-order initialization deliberately match the fixed verification driver;
it is not a second general namelist application. The exact rank and fixed
parity gates establish deterministic configure-time dispatch, not physical
mechanism validity or production MPI scalability.

## Selected normalized mechanism serial reactive EB AMR 2D flow (`0.229.0`)

- [x] configured and installed `pelef_reactive_eb_amr_2d_selected` executable
- [x] shared fixed/selected application driver for static hierarchy setup,
  subcycled advancement, deterministic coarse/fine CSV, volume-weighted
  invariants, and diagnostics
- [x] private selected EB AMR runtime without `pelef_core`, fixed database
  loaders, or either committed generated mechanism
- [x] optional bundle-order composition used by coarse/fine and configured
  physical-boundary initialization
- [x] generated explicit/implicit chemistry policy used by both chemistry
  half-steps on both levels
- [x] independently ordered explicit two-species fixture with 32/8/24 coarse
  and 84/12/48 fine regular/cut/covered cell counts, positive state, closure,
  active H2 change, exact average-down, and repeat-byte identity on both levels
- [x] selected implicit full-H2/O2 coarse and fine outputs byte-identical to
  the fixed executable with active chemistry, all molecular-transport terms,
  characteristic PLM, and an isothermal moving no-slip embedded wall
- [x] unknown species, out-of-range reflected wall temperature, every lexical
  input/coarse/fine alias, dynamic/three-level/multipatch/dynamic-parent modes,
  each checkpoint control, restart request, and fixed-application
  selected-model rejection with input preservation and no level output
- [x] generator namespace reservation plus maximum-length selected-template
  compilation
- [x] eight-configuration 3,478/3,478 regression matrix with Debug, Release,
  Cantera, MPI, and CVODE coverage
- [x] sanitized tests-disabled 567-step Release build/install with 37
  build/install byte-identical ELF executables, non-executable stacks, no
  RPATH/RUNPATH, resolved dependencies, and one coherent OpenMPI ABI
- [x] installed fixed/selected EB AMR Release hash parity plus retained
  fixed/selected regular MPI 1/2/4-rank byte identity
- [ ] selected dynamic/three-level/multipatch EB AMR, checkpoint/restart,
  selected EB 3D, MPI AMR/EB, runtime loading, CFD CVODE, mechanisms above 32
  species, thread safety/performance, detailed fuels, or external validation

The selected target links the selected single-level EB runtime and privately
compiles only the generic hierarchy, reflux, regrid, transport, driver, and
application modules. This evidence closes configure-time dispatch for the
static serial two-level boundary only; it does not qualify a new mechanism
physically or establish external PeleC field parity.

## Selected normalized mechanism serial reactive EB 3D flow (`0.230.0`)

- [x] configured and installable `pelef_reactive_eb_3d_selected` executable
- [x] shared fixed/selected application driver for planar geometry, adaptive
  CFL selection, optional `R-T-H-T-R`, redistribution, diagnostics, and
  dynamic-species CSV output
- [x] private selected EB 3D runtime without `pelef_core`, fixed database
  loaders, or either committed generated mechanism
- [x] optional bundle-order composition used by regular and cut-cell
  initialization and generated explicit/implicit policy used by both masked
  chemistry half-steps
- [x] independently ordered H2/H fixture with 32 regular, 16 cut, and 16
  covered cells, active H2 change, closure, positivity, and repeat-byte
  identity
- [x] fixed wrapper byte-identical to its frozen pre-refactor Debug and
  Release outputs
- [x] complete test-only seven-species selected bundle byte-identical to the
  fixed executable in Debug and Release
- [x] pinned full-H2/O2 selected run with implicit chemistry, viscosity,
  thermal conduction, species diffusion, barodiffusion, 288 regular, 48 cut,
  and 144 covered cells
- [x] invalid, duplicate, nonfinite, zero, and negative compositions;
  out-of-range temperature; input/output alias; checkpoint/restart;
  unsupported element family; and fixed-parser selected-model rejection
- [x] generator dependency/namespace/template compile gates and dynamic
  selected-EB3D result checker
- [x] eight-configuration 3,654/3,654 regression matrix with Debug, Release,
  Cantera, MPI, and CVODE coverage
- [x] sanitized tests-disabled 595-step Release build/install with 38
  build/install byte-identical ELF executables, non-executable stacks, no
  RPATH/RUNPATH, resolved dependencies, and one coherent OpenMPI ABI
- [x] installed selected full-H2/O2 EB 3D checker plus retained fixed EB 3D,
  EB AMR fixed/selected, and regular MPI 1/2/4-rank parity
- [ ] selected checkpoint/restart, arbitrary element diagnostics, general EB
  geometry, EB AMR/MPI 3D, runtime loading, CFD CVODE, mechanisms above 32
  species, thread safety/performance, detailed fuels, or external validation

The selected target links the selected regular-3D runtime and privately
compiles only the generic planar EB geometry, CFL, wall-flux, redistribution,
hydro, transport, driver, and shared application modules. Selected persistence
is rejected before output because the fixed checkpoint schema does not store
the generated bundle fingerprint, composition, or integrator policy. This
evidence qualifies deterministic configure-time dispatch for the bounded
serial planar path; it is not a claim of arbitrary mechanism physics or full
PeleC parity.

## Selected serial reactive EB 3D checkpoint/restart (`0.231.0`)

- [x] byte-preserving fixed schema 1 and exclusive selected schema 2
- [x] selected context records for the configure-time bundle SHA-256,
  generated chemistry-integrator policy, and normalized bundle-order
  composition
- [x] all-or-none selected checkpoint API with no silent schema downgrade
- [x] fixed-reader/schema-2 and selected-reader/schema-1 rejection
- [x] transactional rejection for changed bundle SHA, integrator, composition,
  configuration, mechanism, transport, geometry, state, or diagnostics
- [x] malformed marker/SHA/integrator/count, negative/nonfinite composition,
  truncation, and wrong-shape rollback gates
- [x] selected input/checkpoint and input/restart lexical alias rejection
- [x] separate-process uninterrupted, checkpoint-stop, and restart regression
  with exact 1152-row final CSV and cumulative-diagnostic identity
- [x] selected/fixed uninterrupted reference byte identity in Debug and Release
- [x] automated unchanged fixed Debug/Release schema-1 checkpoint and
  stopped-CSV SHA-256 gates
- [x] 3,686/3,686 eight-configuration matrix and fresh 595-step,
  38-executable install audit with installed fixed/schema-2 restart smoke

The complete selected context is checked before any restart target is
published. The qualified end-to-end case uses the complete seven-species
elementary selected bundle, explicit chemistry, all planar molecular-transport
terms, StateRedist, a first-step stop, and an independent continuation to step
two. A production full-H2/O2 installed smoke also resumes an implicit-policy
schema-2 checkpoint through four hydro steps with exact final output and
cumulative diagnostics. The selected bundle SHA is provenance, not a digest of the checkpoint
body: arbitrary finite, physically consistent payload alteration is not
tamper-evident. General geometry, AMR/MPI EB 3D, distributed or crash-atomic
I/O, schema migration, runtime mechanism loading, CFD CVODE, detailed fuels,
thread safety, and performance remain open.

## Selected static two-level reactive EB AMR 2D restart (`0.232.0`)

- [x] byte-preserving fixed schema 3 and exclusive selected schema 4
- [x] schema-4 context records for bundle SHA-256, generated integrator,
  species count, and normalized bundle-order composition
- [x] all-or-none API with fixed-reader/schema-4 and selected-reader/schema-3
  rejection and no fallback
- [x] transactional rejection for changed bundle SHA, integrator,
  composition, configuration, geometry, wall controls, patch, or fields
- [x] incomplete context, malformed marker, truncation, invalid composition,
  rollback, and non-destructive failed-write unit gates
- [x] actionable checkpoint failure messages propagated through the shared
  application
- [x] input/coarse/fine/checkpoint/restart lexical alias rejection
- [x] separate-process full-H2/O2 fixed reference, selected reference,
  first-step checkpoint-stop, restart, and composition-mismatch regression
- [x] exact fixed/selected and uninterrupted/restarted coarse and fine CSV
  with 64 coarse and 144 fine rows
- [x] frozen selected reference, stopped, and schema-4 checkpoint SHA-256
  gates plus unchanged fixed schema-3 checkpoint and stopped-output hashes
- [x] 3,726/3,726 eight-configuration matrix covering Debug, Release,
  Cantera, MPI, CVODE, and their qualified combinations
- [x] sanitized tests-disabled 595-step Release build/install with 38
  build/install byte-identical ELF executables, non-executable stacks, no
  RPATH/RUNPATH, resolved dependencies, and one coherent OpenMPI ABI
- [x] installed fixed full-H2/O2 reference plus installed selected reference,
  checkpoint-stop, and independent schema-4 restart exact-hash audit
- [ ] selected dynamic, three-level, multipatch, and MPI EB AMR persistence;
  crash-atomic/scalable I/O; payload authentication; runtime loading; CFD
  CVODE; detailed fuels; thread safety; performance; external validation

The qualified selected checkpoint is written after coarse step one at
`6.3479525559440055e-9 s` and resumes to coarse step two at `1.0e-8 s`.
Final selected coarse and fine files match both their uninterrupted selected
references and the fixed references byte-for-byte. Schema 4 is deliberately
limited to the established serial static exactly-two-level hierarchy; it does
not imply persistence support for the other EB AMR topologies.

## Selected dynamic two-level reactive EB AMR 2D restart (`0.233.0`)

- [x] byte-preserving fixed schema 3 and selected static schema 4
- [x] exclusive selected dynamic schema 5 with the schema-4 selected-context
  contract
- [x] finite initial composite-integral baseline with corruption rejection and
  cumulative conservation-diagnostic restart identity
- [x] schema selection bound to fixed/selected and static/dynamic call mode
  with no fallback
- [x] schema-5 round trip plus static/schema-5 and dynamic/schema-4 rejection
  with complete rollback
- [x] selected dynamic front-end support limited to serial, exactly two
  levels, and one fine patch
- [x] continued three-level, multipatch, and dynamic-parent startup rejection
- [x] separate-process full-H2/O2 fixed reference, selected reference,
  topology-changing first-step checkpoint-stop, restart, and composition
  mismatch regression
- [x] exact fixed/selected and uninterrupted/restarted coarse and fine CSV
  with frozen reference, stopped, and schema-5 checkpoint SHA-256 values
- [x] unchanged schema-3/schema-4 restart regression and frozen hashes
- [x] complete eight-configuration Debug/Release, Cantera, MPI, and CVODE
  matrix: 3,766/3,766 tests passed
- [x] sanitized tests-disabled Release build/install/ELF audit with 38/38
  byte-identical ELF executables and installed fixed plus selected dynamic
  restart smoke
- [ ] selected three-level, multipatch, dynamic-parent, and MPI EB AMR
  persistence; crash-atomic/scalable I/O; payload authentication; runtime
  loading; CFD CVODE; detailed fuels; thread safety; performance; external
  validation

The qualified checkpoint follows an actual topology change from coarse patch
`(2:5,2:5)` to `(5:12,4:12)`, records one completed regrid at coarse step one,
and resumes to step three. Schema 5 is deliberately limited to this serial
single-child topology class and does not broaden schema 4.

## Selected static three-level reactive EB AMR 2D restart (`0.234.0`)

- [x] byte-preserving fixed static magic/schema 3 and fixed dynamic
  magic/schema 4
- [x] selected schema 4 under the static magic with complete selected context
  and finite initial composite-integral baseline
- [x] generated bundle-order composition and integrator propagation through
  all root/middle/finest initialization, boundary, and chemistry stages
- [x] transactional reader publication only after context, baseline, topology,
  every state/temperature field, EOS consistency, and terminal marker validate
- [x] fixed/selected cross-schema rejection; bundle, integrator, and
  composition mismatch; corrupt marker/numeric baseline; non-destructive
  invalid write; and empty/default failure-output unit gates
- [x] selected static three-level application/front-end dispatch with continued
  dynamic-three-level, multipatch, and dynamic-parent rejection
- [x] input/root/middle/finest/checkpoint/restart lexical-alias rejection
- [x] separate-process full-H2/O2 fixed reference, selected reference,
  first-step checkpoint-stop, restart, and composition-mismatch chain
- [x] exact fixed/selected and uninterrupted/restarted root, middle, and finest
  CSV with frozen reference, stopped, and checkpoint SHA-256 values
- [x] complete eight-configuration Debug/Release, Cantera, MPI, and CVODE
  matrix: 3,838/3,838 tests passed
- [x] sanitized tests-disabled Release build/install/ELF audit and installed
  fixed plus selected static-three-level restart smoke
- [ ] selected dynamic-three-level, multipatch, dynamic-parent, and MPI EB AMR
  persistence; crash-atomic/scalable I/O; payload authentication; runtime
  loading; CFD CVODE; detailed fuels; thread safety; performance; external
  validation

The selected checkpoint is written after coarse step one at
`6.2673103087327447e-9 s` and resumes to coarse step two at `1.0e-8 s`.
Selected root/middle/finest final files match the uninterrupted selected and
fixed references byte-for-byte, and restart reproduces the uninterrupted
`4.5441946758675300e-12` maximum composite conservation error text. This
schema is limited to the serial static one-patch-per-level topology.

## Selected dynamic three-level reactive EB AMR 2D restart (`0.235.0`)

- [x] byte-preserving fixed dynamic-three-level magic/schema 4
- [x] selected schema 5 under the dynamic magic with complete selected context
  and finite original composite-integral baseline
- [x] selected application/front-end support for the serial exactly-three-
  level, one-patch-per-level dynamic lifecycle, including dynamic-parent
  regridding
- [x] generated bundle-order composition and integrator propagation through
  every root/middle/finest initialization, boundary, and chemistry stage
- [x] transactional publication only after context, baseline, dynamic policy,
  committed patches, every field, EOS consistency, terminal marker, and end-
  of-stream validation
- [x] fixed/selected cross-schema and static/dynamic cross-magic rejection with
  empty/default failure targets
- [x] baseline closure/marker/numeric corruption, dynamic controls, middle
  patch, first field, and terminal-marker corruption gates plus non-destructive
  invalid writes
- [x] repurposed negative configuration gates for finest-patch removal and
  dynamic-parent without dynamic-three-level regridding
- [x] separate-process full-H2/O2 fixed reference, selected reference,
  dynamic-parent topology-changing first-step checkpoint-stop, independent
  restart, and composition-mismatch chain
- [x] exact fixed/selected and uninterrupted/restarted root, middle, and finest
  CSV with frozen reference, stopped, and schema-5 checkpoint SHA-256 values
- [x] complete eight-configuration Debug/Release, Cantera, MPI, and CVODE
  matrix: 3,886/3,886 tests passed
- [x] sanitized tests-disabled Release build/install/ELF audit and installed
  fixed plus selected dynamic-three-level restart smoke
- [ ] selected multipatch and MPI EB AMR persistence; crash-atomic/scalable
  I/O; payload authentication; runtime loading; CFD CVODE; detailed fuels;
  thread safety; performance; external validation

The qualified lifecycle moves the middle patch from `(2:11,2:11)` to
`(5:10,3:10)` and the finest patch from `(6:9,6:9)` to `(3:10,3:14)` before
its first-step checkpoint. It resumes at coarse step one and reaches step two
at `1.0e-8 s`. Selected root/middle/finest final files match both the
uninterrupted selected and fixed references byte-for-byte, and restart
reproduces the uninterrupted `1.1879823316399355e-10` maximum composite
conservation error text. This schema is limited to the serial dynamic
one-patch-per-level topology.

## Selected dynamic multipatch reactive EB AMR 2D restart (`0.236.0`)

- [x] byte-preserving fixed multipatch magic/schema 3, pinned by checkpoint
  SHA-256 `3ce931f5dcd017de4864789f53b57a47e40fe82cc5fa6c87e9292ff7efe1022c`
- [x] selected schema 4 under the same magic with complete selected context
  and finite original composite-integral baseline
- [x] selected application/front-end support for the serial dynamic exactly-
  two-level multipatch lifecycle
- [x] generated bundle-order composition and integrator propagation through
  root and every child initialization, boundary, and chemistry stage
- [x] transactional publication only after context, baseline, committed patch
  set, every field, EOS consistency, terminal marker, and end-of-stream
  validation
- [x] fixed/selected cross-schema rejection; composition mismatch; corrupt
  baseline, terminal marker, and trailing-content rejection; non-destructive
  invalid writes; and empty/default failed-read targets
- [x] selected front-end lexical-alias rejection for every possible derived
  `_patchNNNN` child output before mechanism loading or file creation
- [x] separate-process full-H2/O2 fixed reference, selected reference,
  first-step checkpoint-stop, independent restart, and composition-mismatch
  chain with two fine patches
- [x] exact fixed/selected and uninterrupted/restarted root, patch 1, and
  patch 2 CSV with frozen reference, stopped, and schema-4 checkpoint SHA-256
  values
- [x] complete eight-configuration Debug/Release, Cantera, MPI, and CVODE
  matrix: 3,934/3,934 tests passed
- [x] sanitized tests-disabled Release build/install/ELF audit with 38/38
  byte-identical ELF executables and installed fixed plus selected multipatch
  restart smoke
- [ ] selected MPI EB AMR persistence; crash-atomic/scalable I/O; payload
  authentication; runtime loading; CFD CVODE; detailed fuels; thread safety;
  performance; external validation

The qualified serial case uses a 14-by-14 root and two disjoint 10-by-10
children with coarse-cell bounds `(2:6,6:10)` and `(9:13,9:13)`. Its selected
checkpoint is written after coarse step one at
`4.5125253966820509e-9 s` and resumes to coarse step three at `1.0e-8 s`.
Selected root and both child files match the uninterrupted selected and fixed
references byte-for-byte; the maximum mass-fraction closure error is
`2.2204460492503131e-16`. This schema is limited to the established serial
dynamic two-level patch-set topology and does not imply MPI persistence.

## Selected sparse MPI reactive AMR 1D restart (`0.237.0`)

- [x] byte-preserving fixed patch-tree magic/schema 1 with pinned checkpoint
  SHA-256 `1337b54d2761f2e3a2d25e264e9d6d073d1934756b206da16634d3140d933ed6`
- [x] selected schema 2 with bundle SHA-256, generated integrator,
  bundle-order composition, and original composite-integral baseline
- [x] installed configure-time selected sparse MPI AMR 1D front end and shared
  fixed/selected MPI application lifecycle
- [x] communicator-wide context consensus before initialization or output
- [x] generated integrator propagation through every sparse chemistry half-
  step, with collective invalid-policy rollback on one, two, and four ranks
- [x] transactional schema/context/baseline/payload/end-of-stream validation;
  cross-schema, truncation, trailing-content, and mismatch rejection
- [x] communicator-neutral checkpoint without rank count or owner map, with
  one-rank stop and independent two-/four-rank restart
- [x] exact fixed/selected, one-/two-/four-rank, and uninterrupted/restarted
  final CSV bytes for the pinned full-H2/O2 three-level process
- [x] bundle, integrator, composition, species-order, baseline, rank-context,
  and input/output-alias failure gates without output publication
- [x] density, raw species, near-floor negative species, closure, temperature,
  and geometry corruption rejection before selected restart publication
- [x] final-source MPI Debug focused set: 37/37 tests passed; independent Luna
  Max read-only audit reports no blocking defect
- [x] complete eight-configuration Debug/Release, Cantera, MPI, and CVODE
  matrix: 4,020/4,020 tests passed
- [x] sanitized tests-disabled Release build/install/ELF audit with 39/39
  byte-identical ELF executables and installed fixed plus selected one-/two-/
  four-rank reference, changed-rank restart, and 12-case failure smoke
- [ ] selected MPI EB AMR; crash-atomic/scalable I/O; payload authentication;
  runtime loading; CFD CVODE; detailed fuels; thread safety; performance;
  external validation

The selected checkpoint is written after coarse step one at
`2.8740562104126945e-10 s` and resumes to coarse step four at
`1.0000000000000001e-9 s`. All final CSV files have SHA-256
`664fd578e102ddb769972289bd60e3c2678b61b5d4b9b219d053a682a90ca809`;
the schema-2 checkpoint SHA-256 is
`0e5cc002ac5cceb3fa7732f8fe46382f140e038d19cac17664230ce8c500c2fd`.
This is a sparse MPI AMR 1D contract and does not imply selected MPI EB AMR.

## Selected sparse MPI reactive EB patch-tree 2D (`0.238.0`)

- [x] shared fixed/selected mechanism-independent MPI EB patch-tree driver
- [x] installed configure-time selected MPI reactive EB patch-tree executable
- [x] complete bundle SHA-256, normalized composition, species/reaction/
  transport, and generated-integrator context validation
- [x] communicator-wide selected-context consensus before root initialization
- [x] owner-only selected root initialization and selected boundary composition
- [x] generated integrator propagation through public chemistry, full-physics,
  and to-time entry points and both chemistry half-steps
- [x] collective invalid-policy and rank-policy rejection with exact sparse
  state non-mutation on one, two, and four ranks
- [x] four-level solution-driven dynamic regridding with active full-H2/O2
  chemistry and exact fixed-one-rank/selected-one-/two-/four-rank output
- [x] two-species fixture-bundle execution and CSV validity gate for the
  generated explicit-integrator path
- [x] selected model, fixed model, input/output alias, checkpoint cadence,
  checkpoint path, and restart path rejection without output publication
- [x] fixed schema-8 GNU Debug/Release checkpoint hash gates and Release
  byte-for-byte comparison with the prior `0.237.0` executable
- [ ] selected MPI EB checkpoint/restart schema; three-level or multipatch
  fixed-depth MPI execution; runtime loading; scalable/crash-atomic I/O;
  detailed fuels; performance; physical or external PeleC validation

The public parity case advances a 12-by-12 planar-EB root to `1.0e-12 s`,
performs dynamic regridding at initialization and after its first coarse step,
and activates levels 0 through 3. Fixed one rank and selected one, two, and
four ranks produce the same composite CSV SHA-256
`96cfbb523031b8a22163ff8308627ba937a7d3789719de2c81011e50490009dc`.
The frozen GNU Release schema-8 checkpoint SHA-256 is
`d3ed831ad70a488f7303bbc296b07f63f503c9f02c04a0ddff81d0bb8e7f1de1`.

## Selected sparse MPI reactive EB patch-tree 2D restart (`0.239.0`)

- [x] fixed magic/schema 8 byte preservation and fixed/selected cross-schema
  rejection
- [x] exclusive selected schema 9 with bundle SHA-256, generated integrator,
  bundle-order composition, complete fingerprint, and original composite
  baseline
- [x] selected write validation before target open and private transactional
  read validation through `END_CHECKPOINT` plus strict end-of-stream
- [x] finite positive density/energy/temperature, nonnegative raw species,
  species-closure, EOS-recovery, geometry, clock, counter, and regrid-history
  validation
- [x] communicator-wide exact selected-context consensus in the MPI I/O API
- [x] collective selected metadata-presence rejection before gather/root read:
  all eight writer arguments and the reader fingerprint, including direct
  two-rank read/write mismatch and no-publication coverage
- [x] rank-neutral topology without communicator size or owner map, with
  one-rank stop and independent two-/four-rank restart
- [x] exact fixed/selected and uninterrupted/restarted final CSV bytes on one,
  two, and four ranks
- [x] schema, species order, bundle, integrator, composition, fingerprint,
  geometry, baseline, density, raw/near-floor species, closure, temperature,
  end-marker, truncation, and trailing-content rejection without output
- [x] all six lexical input/output/checkpoint/restart alias gates before
  mechanism loading or file publication
- [x] GNU Debug and Release schema-9 checkpoint/final/stopped SHA-256 gates
- [x] final full configuration matrix: 4074/4074 tests passed; configured
  full-H2/O2 selected MPI Debug tree: 573/573 tests passed
- [x] tests-disabled clean Release install: 40/40 ELF audit plus installed
  fixed/selected one-/two-/four-rank restart, 16 corruption, and six alias gates
- [ ] crash-atomic/scalable I/O; payload authentication; schema migration;
  fixed-depth MPI modes; runtime loading; CFD CVODE; detailed fuels;
  performance; physical or external PeleC validation

The public restart case uses an 8-by-8 planar-EB root and four dynamically
active levels. It writes schema 9 after coarse step one at
`5.1166209105049765e-13 s` and resumes to coarse step two at `1.0e-12 s`.
The schema-9 checkpoint SHA-256 is
`2edd758ca1ff888735f5719de8ed9969b36f60d45bb0b3c523e47c2e35510f4d`;
all fixed, selected, and restarted final CSV files have SHA-256
`5394fd007283d746ec4507832ccce8a7013fd4fb9b5c1b10d50d7c91e216975a`.

## Rank-local sparse static 3D AMR hydro (`0.240.0`)

- [x] uneven rank-local coarse x slabs and parent-aligned fine ownership
- [x] valid empty fine payloads and active previous/next fine-rank discovery
- [x] hierarchy storage limited to local conserved state and temperature
- [x] two-plane periodic coarse halo propagation for one-cell slabs
- [x] active-neighbor fine halo exchange plus time-interpolated coarse/fine
  exterior ghosts
- [x] rank-local PCM and characteristic-PLM SSPRK2 on coarse and fine levels
- [x] deterministic six-face fine flux accumulation, local reflux,
  average-down, and temperature recovery
- [x] transactional collective commit and rank-dependent metadata rejection
- [x] ownership, patch, NASA7, floating-control, solver, reconstruction, and
  limiter fingerprints before data-dependent collectives
- [x] root-only rank-neutral scatter/gather retaining checkpoint schema 2
- [x] exact serial versus MPI coarse/fine output at one, two, four, and eight
  ranks for PCM and characteristic PLM
- [x] two-rank checkpoint-stop with independent four-/eight-rank exact restart
- [x] direct coverage for two-plane halo bits, zero-fine ranks, fine-neighbor
  ranks, sparse storage reduction, rollback, rank-dependent layouts, and
  collective rejection of invalid coarse state on a zero-fine rank
- [x] complete eight-configuration Debug/Release, Cantera, MPI, and CVODE
  matrix: 4096/4096 tests passed
- [x] tests-disabled clean Release install: 40/40 ELF byte identity, dependency
  audit, installed 1/2/4/8-rank PLM parity, and changed-rank restart parity
- [ ] dynamic topology; AMR chemistry/transport; physical coarse boundaries;
  CTU/PPM; 3D EB AMR; scalable I/O; performance and external field parity

The public entropy-wave outputs retain serial hashes
`864e9cabe8f5c2bee19f97e9279bdada7647305171955c6c767cc83e4db15c21`
for the coarse level and
`8f764e38491e9cf62f8c0f17eac2470c4705f7d80bfe7e686db4ba54de2f58d7`
for the fine level under characteristic PLM. The checkpoint stores neither
rank count nor ownership. Formatted root gather/scatter is compatibility I/O
and is not a distributed-I/O claim.

## Fixed chemistry in static 3D AMR (`0.241.0`)

- [x] serial elementary/full-H2/O2 coarse and fine cell-local source phases
- [x] source-phase average-down and covered-coarse temperature recovery
- [x] hierarchy-transactional `R(dt/2)-H(dt)-R(dt/2)` composition
- [x] exact chemistry-disabled dispatch through the `0.240.0` hydro path
- [x] rank-local sparse MPI source integration with valid zero-fine ranks
- [x] communicator-wide ordered reaction, tolerance, integrator, solver,
  reconstruction, limiter, chemistry-enable, NASA7, patch, and ownership
  consensus
- [x] collective rollback for invalid integrators, reaction mismatch,
  rank-dependent chemistry policy, and a failed coarse source on a zero-fine
  rank
- [x] public full-H2/O2 serial control/reacting cases with strict CSV schema,
  topology, physicality, species closure, activity, Euler, and H/O/N gates
- [x] exact serial versus MPI public coarse/fine bytes at one, two, four, and
  eight ranks
- [x] chemistry-enabled schema-2 checkpoint/restart rejection before execution
- [x] complete eight-configuration Debug/Release, Cantera, MPI, and CVODE
  matrix: 4166/4166 tests passed
- [x] tests-disabled clean Release install: 40/40 ELF byte identity and
  dependency audit, installed serial/MPI 1/2/4/8 chemistry parity, and
  changed-rank characteristic-PLM restart parity
- [x] independent Luna Max final read-only audit: zero release-blocking issues
- [ ] configure-time selected mechanisms; molecular transport; dynamic 3D
  topology; physical coarse boundaries; CTU/PPM; 3D EB AMR; distributed I/O;
  performance and external PeleC field parity

The full-H2/O2 hotspot at `1.0e-9 s` has checker species-change metric
`6.481e-6` and changes the fine temperature by `2.066e-4 K` relative to its
inert control. The checker reports Euler error `1.320e-16` and H/O/N error
`4.377e-16`. Reacting serial and one-/two-/four-/eight-rank MPI files have
SHA-256 `3313a9a6eff5d6175625e4bbe3613a782cdd794314fd9fbdd53ddaa309fc9516`
for coarse and
`4ff0785c3462a49648fd694d090eb23bffd5af05a43a33360cd187417bd8a31d`
for fine output. This is a short structural and deterministic parity case,
not ignition, detailed-fuel, performance, or physical validation evidence.

## Selected mechanisms in static 3D AMR (`0.242.0`)

- [x] explicit `thermo_model='selected'` opt-in for the serial and sparse-MPI
  static two-level 3D AMR applications
- [x] mechanism-independent serial and MPI lifecycle drivers shared by the
  fixed and generated frontends
- [x] exact species-name composition, bundle SHA-256, ordered NASA7,
  reaction, gas-transport, and generated integrator provenance
- [x] pre-initialization communicator consensus with representation-exact
  floating-point comparison and collective rejection without published state
- [x] zero-fine-rank participation in context, source, and acceptance
  collectives
- [x] generated two-species fixture parity in serial and at one, two, and four
  MPI ranks, including a four-rank zero-fine ownership case
- [x] generated full-H2/O2 parity with the fixed serial and one-, two-, four-,
  and eight-rank MPI applications
- [x] startup rejection of selected transport, checkpoint, restart, and
  lexical input/output alias requests; the fixed frontend still rejects the
  selected model
- [x] fixed full-H2/O2 and characteristic-PLM public artifact SHA-256 gates
- [x] complete eight-configuration matrix: `4320/4320` tests passed
- [x] clean tests-disabled Release: `659/659` build steps, 42/42 installed ELF
  audit, selected serial/MPI 1/2/4/8 parity, and fixed changed-rank PLM restart
- [ ] selected context-bound static 3D AMR restart; AMR molecular transport;
  dynamic 3D topology; physical coarse boundaries; CTU/PPM; 3D EB AMR;
  distributed I/O; performance and external PeleC field parity

This increment qualifies deterministic configure-time dispatch and exact
serial/rank parity. It does not qualify ignition, detailed-fuel physics,
production scaling, or external PeleC field parity.

## Selected static 3D AMR restart (`0.243.0`)

- [x] fixed schema-2 checkpoint bytes and changed-rank restart gates preserved
- [x] exclusive schema 3 requiring complete selected context before file open
- [x] bundle SHA-256, generated integrator, normalized composition, complete
  ordered reaction records, and chemistry tolerances serialized and validated
- [x] ordered NASA7, numerical/topology fingerprint, original baseline,
  clock/reflux history, both levels, terminal marker, and strict EOF retained
- [x] fixed/selected cross-schema rejection and transactional context mismatch
  rejection without target-state publication
- [x] communicator-wide context consensus hardened against malformed
  rank-local third-body-efficiency shapes before serialization
- [x] chemistry tolerances included in exact context consensus before root I/O
- [x] lexical input/output/checkpoint/restart collision gates before file write
- [x] serial full-H2/O2 uninterrupted/checkpoint/restart exact parity
- [x] two-rank checkpoint resumed exactly on one, two, four, and eight ranks,
  including zero-fine-plane ranks
- [x] byte-identical serial/MPI checkpoint and frozen schema-3/final hashes
- [x] final eight-configuration matrix: `4366/4366` tests passed
- [x] clean tests-disabled Release: `659/659` build steps, 42/42 installed ELF
  audit, selected serial checkpoint resumed exactly in serial and at MPI
  ranks 1/2/4/8, and fixed schema-2 MPI 2-to-4/8 restart parity
- [x] independent Luna Max final audit: zero release-blocking findings
- [ ] payload authentication; crash-atomic/scalable I/O; schema migration; AMR
  transport; dynamic topology; physical boundaries; 3D EB AMR; performance;
  detailed-fuel physical or external PeleC validation

The selected schema-3 checkpoint SHA-256 is
`f44f69a7e1a1657b549b26373feb277840732e5453b0734321cc8c7cdf91c423`.
Exact restarted coarse/fine CSV SHA-256 values are
`a03884837909bab719f0783d23ea9150a81f2afc628fee1a01c645c1b47e0ff3` and
`8d5a992d07b89832624e37af6ddacd510df63633707678693b9053057d4efcbe`.

## Molecular transport in static 3D AMR (`0.244.0`)

- [x] coarse/fine parabolic stable-step reduction with `r^2` fine scaling
- [x] serial coarse SSPRK2 and `r^2`-subcycled fine SSPRK2 transport
- [x] time-interpolated, limited-PLM transport ghosts on all patch faces
- [x] viscosity, thermal conduction, mixture-averaged species diffusion,
  barodiffusion, correction velocity, and species-enthalpy flux
- [x] time- and area-averaged x/y/z lower/upper interface flux registers
- [x] six-face uncovered-coarse species limiter with a conservative fine-side
  flux correction, reflux, average-down, and EOS temperature recovery
- [x] hierarchy-transactional `R-T-H-T-R` with late hydro rollback
- [x] rank-local sparse MPI transport with valid zero-fine-plane ranks
- [x] exact transport metadata and policy consensus before payload collectives
- [x] direct invalid external-theta and barodiffusion-dependency rejection
- [x] serial adversarial limiter activation on all six faces, individual
  process modes, refinement ratio three, and composite conservation
- [x] fixed and selected public full-H2/O2 serial transport identity
- [x] fixed and selected one-/two-/four-/eight-rank MPI exact parity with
  frozen coarse/fine SHA-256 values
- [x] explicit transport checkpoint/restart rejection before output
- [ ] dynamic topology; boundary-touching refinement; physical coarse
  boundaries; transport persistence; 3D EB AMR; scalable I/O; performance;
  detailed-fuel physical or external PeleC validation

The final eight-configuration qualification passes `4432/4432` tests. The
tests-disabled clean Release installs and audits `42/42` ELF executables, then
runs installed fixed/selected transport in serial and at one, two, four, and
eight ranks. Installed selected schema-3 and fixed schema-2
transport-disabled restart compatibility also remains exact. See
[`validation/0.244.0.md`](validation/0.244.0.md) and its frozen evidence.

The qualified public case ends at `1.0e-9 s`; it is a short structural,
conservation, and deterministic-parity gate. It does not establish ignition,
detailed-fuel accuracy, production scaling, or whole-field PeleC parity.

## Static 3D AMR transport checkpoint persistence (`0.245.0`)

- [x] exclusive schema 4 for fixed public `full_h2o2` transport with
  `chemistry_enabled=.false.` and `transport_enabled=.true.`
- [x] exclusive schema 5 for selected chemistry with
  `thermo_model='selected'`, `chemistry_enabled=.true.`, and
  `transport_enabled=.true.`
- [x] complete ordered gas-transport database, parameter convention, operator
  identity, transport controls, and cumulative transport diagnostics bound to
  transport-enabled checkpoints
- [x] `POST_ACCEPTED_COARSE_STEP` checkpoint boundary after committed
  `R(dt/2)-T(dt/2)-H(dt)-T(dt/2)-R(dt/2)` transport/chemistry/hydro stepping
- [x] private candidate read, strict terminal marker/EOF and close validation,
  and publication only after every checkpoint check succeeds
- [x] rank-neutral root-formatted serial/MPI checkpoint with all-rank context
  consensus, zero-fine-plane participation, collective root-I/O failure, and
  exact changed-rank restart
- [x] fixed schema-2 and selected schema-3 transport-disabled compatibility
  retained without cross-schema fallback
- [x] focused GNU Debug/Release serial transport restart gates: `10/10` each
- [x] focused GNU/OpenMPI Debug/Release MPI transport restart gates: `23/23`
  each
- [ ] full configuration matrix, clean installed Release qualification,
  crash-atomic replacement, payload authentication, scalable I/O, dynamic
  topology, physical boundaries, 3D EB AMR, performance, and external PeleC
  field parity

The focused schema-4/schema-5 checkpoint, coarse-output, and fine-output
SHA-256 values are recorded in
[`validation/0.245.0.md`](validation/0.245.0.md). No 0.245 full-matrix or
clean-install result is claimed; the `0.244.0` counts and install evidence
above remain historical and are not reused for this increment.
