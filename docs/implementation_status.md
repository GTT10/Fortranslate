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
- [x] qualified general-EOS PeleC-style Riemann parity
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
