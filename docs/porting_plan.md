# PeleC → Modern Fortran Porting Plan

## Project goal

Reimplement the numerical algorithms, physics capabilities, and regression behavior of PeleC in Modern Fortran as an independent solver. This is not intended to be a mechanical C++→Fortran syntax translation. The implementation should be developed as a clean, testable Fortran codebase with explicit parity checks against PeleC.

Working name: **PeleF**.

Target end state:

- compressible reacting-flow solver
- multispecies transport and chemistry
- MPI parallelism
- AMR
- embedded boundaries
- LES
- Lagrangian particles / spray
- restart and visualization output
- CPU and eventually GPU execution

PeleC reference branch: `Pele-Suite/PeleC:development`.

Development branch in this repository: `agent/pelec-fortran-port`.

---

## Core principles

1. Port behavior and algorithms, not C++ syntax.
2. Do not attempt to recreate all of AMReX at once.
3. Implement only the infrastructure required by the PeleC feature currently being ported.
4. Every major capability must have an automated parity gate before the next layer is added.
5. Keep the solver usable at intermediate stages.
6. Separate physics kernels from infrastructure and parallel execution.
7. Prefer simple Modern Fortran modules and explicit array kernels over heavy OOP.
8. Preserve numerical traceability: every Fortran subsystem should map to a corresponding PeleC subsystem and regression case.

---

## Repository documentation

Create and maintain:

```text
docs/
├── architecture.md
├── pelec_mapping.md
├── dependency_map.md
├── state_variables.md
├── numerical_methods.md
├── parity_strategy.md
├── implementation_status.md
└── design_decisions/
```

`pelec_mapping.md` should map reference files and responsibilities, for example:

| PeleC | PeleF |
|---|---|
| `Source/main.cpp` | `app/main.F90` |
| `Source/Advance.cpp` | `src/driver/time_integrator_mod.F90` |
| `Source/Hydro.cpp` | `src/hydro/hydro_mod.F90` |
| `Source/Godunov.*` | `src/hydro/godunov_mod.F90` |
| `Source/Riemann.H` | `src/hydro/riemann_mod.F90` |
| `Source/PLM.H` | `src/hydro/reconstruction_plm_mod.F90` |
| `Source/PPM.*` | `src/hydro/reconstruction_ppm_mod.F90` |
| `Source/WENO.H` | `src/hydro/reconstruction_weno_mod.F90` |
| `Source/Diffusion.*` | `src/diffusion/` |
| `Source/React.cpp` | `src/chemistry/chemistry_integrator_mod.F90` |
| `Source/LES.*` | `src/les/` |
| `Source/PeleCAmr.*` | `src/amr/` |
| `Source/EB.*` | `src/eb/` |
| `Source/Particle.cpp` | `src/particles/` |

---

## Proposed Fortran architecture

```text
app/
└── pelef.F90

src/
├── core/
│   ├── precision_mod.F90
│   ├── constants_mod.F90
│   ├── state_indices_mod.F90
│   ├── mesh_mod.F90
│   └── field_mod.F90
├── physics/
│   ├── eos_mod.F90
│   ├── thermo_mod.F90
│   ├── transport_mod.F90
│   └── chemistry_mod.F90
├── hydro/
│   ├── primitive_mod.F90
│   ├── riemann_mod.F90
│   ├── godunov_mod.F90
│   ├── reconstruction_plm_mod.F90
│   ├── reconstruction_ppm_mod.F90
│   └── reconstruction_weno_mod.F90
├── diffusion/
├── chemistry/
├── les/
├── amr/
├── eb/
├── particles/
├── parallel/
├── io/
└── driver/
```

Implementation style:

- standard Fortran modules
- explicit interfaces
- contiguous arrays
- ordinary `do` loops or `do concurrent`
- no unnecessary derived-type hierarchy
- design kernels so they can later be used with OpenMP/OpenACC/OpenMP target

---

## State vector design

PeleC-style conserved quantities should be represented explicitly, with Fortran-native 1-based indices.

Initial base state layout:

```fortran
integer, parameter :: IRHO = 1
integer, parameter :: IMX  = 2
integer, parameter :: IMY  = 3
integer, parameter :: IMZ  = 4
integer, parameter :: IET  = 5
integer, parameter :: IEI  = 6
integer, parameter :: ITEM = 7
```

Species, advected scalars, auxiliary variables, and future soot variables are appended after the thermodynamic state.

A documented mapping from PeleC indices to PeleF indices must be maintained and tested.

---

# Development phases

## Phase 0 — project skeleton and parity infrastructure

Deliverables:

- CMake build
- gfortran support
- basic test harness
- GitHub Actions
- documentation skeleton
- reference-data layout
- comparison tooling

No CFD functionality required yet.

Parity infrastructure should calculate, where applicable:

- L1 error
- L2 error
- L∞ error
- global mass error
- momentum error
- energy error
- species conservation error

---

## Phase 1 — minimal compressible Euler solver

Scope:

- uniform Cartesian mesh
- serial execution
- no AMR
- no chemistry
- no viscosity
- single ideal gas

Implement:

- conserved ↔ primitive conversion
- constant-γ ideal-gas EOS
- CFL timestep
- Rusanov / local Lax-Friedrichs flux
- finite-volume update
- boundary conditions
- RK2 / method-of-lines stepping

Initial regression cases:

- Sod shock tube
- Shu–Osher
- Sedov blast
- smooth advection
- isentropic vortex

Acceptance gate:

- builds with gfortran
- runs Sod to completion
- positive density and pressure
- bounded conservation errors
- automated analytical/reference comparison

---

## Phase 2 — PeleC-style Godunov hydrodynamics

Implement incrementally:

1. piecewise constant reconstruction
2. PLM
3. PPM
4. characteristic tracing
5. approximate PeleC-style Riemann solver
6. transverse corrections
7. multidimensional Godunov update
8. WENO

Do not implement all methods at once. Every reconstruction level must re-run the same hydro regression suite.

The `0.68.0` milestone extends item 5 from the constant-`gamma` reduction to
the NASA7 ideal-gas-mixture reactive state. It retains the upstream acoustic
star estimate, species-density correction, wave interpolation, and EOS energy
reconstruction behind a separately selectable `pelec` path.

Primary parity cases:

- Sod
- Shu–Osher
- Sedov

---

## Phase 3 — multispecies advection

Add conserved species densities `rho*Y_k`.

Implement:

- species advection
- mass-fraction conversion
- species positivity handling
- normalization / consistency
- `rho = sum(rhoY_k)` consistency checks

Primary parity case:

- PeleC `MultiSpecSod`

At the end of this phase the solver is multispecies but non-reacting and inviscid.

---

## Phase 4 — thermodynamics and EOS

Replace the constant-γ model with a reusable thermodynamics layer.

Order:

1. constant-γ ideal gas
2. NASA-7 polynomial thermodynamics
3. multispecies ideal-gas EOS
4. internal-energy → temperature inversion
5. mixture cp/cv/gamma/sound speed
6. optional SRK EOS later

Required API concepts include:

- `(rho, p, Y) -> c`
- `(rho, p, Y) -> e`
- `(rho, T, Y) -> p`
- `(rho, e, Y) -> T`

---

## Phase 5 — chemistry mechanism code generation

Do not manually rewrite large PelePhysics chemistry mechanisms.

Build a generator:

```text
CHEMKIN / YAML
      ↓
parser
      ↓
internal mechanism representation
      ↓
Fortran generator
      ↓
mechanism_mod.F90
```

Generated data/functions should include:

- molecular weights
- NASA coefficients
- stoichiometry
- Arrhenius rates
- reversible rates
- third-body reactions
- falloff / Troe
- production rates
- optional analytic/sparse Jacobian generation

Mechanisms should be replaceable without changing solver source.

---

## Phase 6 — zero-dimensional reactor

Before coupling chemistry to CFD, implement standalone reactors:

- constant-volume reactor
- constant-pressure reactor

Compare:

- PelePhysics/PeleC
- Cantera
- PeleF

Metrics:

- temperature history
- species histories
- heat-release history
- ignition delay

Primary PeleC regression reference:

- `zeroD`

No 3D chemistry work proceeds until this phase passes.

---

## Phase 7 — stiff chemistry integration

Expose chemistry through a stable API such as:

```fortran
call reactor_integrate(state, dt, reaction_source)
```

Initial integrator strategy:

- use SUNDIALS through its official Fortran interface
- keep the solver independent of PeleC, AMReX, and PelePhysics

A native alternative integrator may be added later, but is not required for the initial port.

---

## Phase 8 — hydro / diffusion / reaction coupling

Keep the subsystems separate and add coupling in increasing complexity.

Order:

1. Strang splitting
2. MOL RK2
3. PeleC-compatible coupled MOL behavior
4. SDC later

Subsystems:

- hydro
- diffusion
- reactions
- external/source terms

---

## Phase 9 — diffusion / Navier–Stokes

Implement:

- viscosity
- thermal conductivity
- species diffusion
- enthalpy diffusion
- prescribed zero-net-mass wall species flux as the transport contract for
  later catalytic surface kinetics (`0.69.0`)

Transport properties must be isolated behind a transport API.

Primary regression cases:

- Taylor–Green vortex
- ChannelFlow
- PMF

---

## Phase 10 — reacting-flow integration

Combine:

- compressible flow
- multispecies thermodynamics
- transport
- detailed chemistry

Primary regression cases:

- zeroD
- TGReact
- PMF

Completion of this phase defines a useful **uniform-grid CPU reacting-flow solver**.

---

## Phase 11 — MPI domain decomposition

Use the standard Fortran MPI interface:

```fortran
use mpi_f08
```

Initial strategy:

- Cartesian block decomposition
- local patches
- halo exchange
- global reductions

Required parity checks:

- 1 rank
- 2 ranks
- 4 ranks
- 8 ranks

Parallel decomposition must not materially change the numerical result.

The `0.20.0` through `0.24.0` milestones complete this initial one-dimensional
slice: uneven block ownership, periodic halos, global reductions, ordered
gather, conservative multispecies hydro, adaptive implicit full-H2/O2
chemistry, molecular transport, and transactional reactive Strang splitting.
The Debug and Release gates exercise 1, 2, 4, and 8 ranks. Multidimensional
decomposition and load balancing are deferred until the AMR data model exists.

---

## Phase 12 — AMR

Do not reimplement all of AMReX. Implement only the AMR subset required by PeleF.

Required concepts:

- level
- patch / box
- refinement ratio
- ghost cells
- prolongation
- restriction
- subcycling
- coarse/fine synchronization
- flux register
- reflux
- tagging
- regridding

Implementation order:

1. [x] static two-level hierarchy (`0.25.0`)
2. [x] restriction (`0.25.0`)
3. [x] prolongation (`0.25.0`)
4. [x] level subcycling (`0.25.0`)
5. [x] reflux (`0.25.0`)
6. [x] tagging (`0.26.0`)
7. [x] dynamic regrid (`0.26.0`)
8. [x] arbitrary multiple levels in hierarchy primitives (`0.30.0`)

Reactive application integration:

- [x] two-level reactive state ownership (`0.27.0`)
- [x] fine subcycling and coarse-time ghost interpolation (`0.27.0`)
- [x] reactive flux-register synchronization (`0.27.0`)
- [x] hierarchy-wide chemistry splitting and rollback (`0.27.0`)
- [x] limited PLM coarse/fine reconstruction (`0.28.0`)
- [x] characteristic PPM coarse/fine reconstruction (`0.34.0`)
- [x] hybrid WENO5-JS/WENO5-Z coarse/fine reconstruction (`0.35.0`)
- [x] WENO3-Z/WENO7-Z coarse/fine reconstruction (`0.36.0`)
- [x] AMR molecular transport (`0.29.0`)
- [x] arbitrary-depth reactive state ownership and recursive advancement
  (`0.31.0`)
- [x] dynamic arbitrary-depth reactive regridding and composite output
  (`0.32.0`)
- [x] overlap-preserving transfer for changed multilevel hierarchies (`0.33.0`)
- [x] outflow physical-boundary refinement and one-sided reflux (`0.37.0`)
- [x] two-level multipatch geometry, regrid transfer, and reactive hydro
  (`0.38.0`)
- [x] two-level multipatch chemistry and molecular transport (`0.39.0`)
- [x] tag-driven two-level multipatch application, regridding, and output
  (`0.40.0`)
- [x] arbitrary-depth multipatch tree geometry and conservative field
  operations (`0.41.0`)
- [x] arbitrary-depth parent-owned flux registers and tree synchronization
  (`0.42.0`)
- [x] arbitrary-depth branched reactive PCM hydro and recursive subcycling
  (`0.43.0`)
- [x] arbitrary-depth patch-tree chemistry splitting and rollback (`0.44.0`)
- [x] arbitrary-depth patch-tree molecular transport and full physics split
  (`0.45.0`)
- [x] plan-driven dynamic patch-tree rebuild and overlap transfer (`0.46.0`)
- [x] per-parent automatic tag clustering and arbitrary-depth patch-tree
  rebuild (`0.47.0`)
- [x] adjacent patch-tree sibling exchange, shared fine/fine flux ownership,
  and reflux suppression (`0.48.0`)
- [x] deterministic MPI patch owner maps, collective hierarchy consensus,
  owner-authoritative field synchronization, and cross-rank sibling halos
  (`0.49.0`)
- [x] owner-only patch-tree chemistry, global acceptance/rollback, and serial
  reacting-field parity (`0.50.0`)
- [x] owner-only recursive patch-tree hydro, shared fine/fine flux handling,
  reflux, rollback, and serial PCM/PPM parity (`0.51.0`)
- [x] owner-only recursive patch-tree molecular transport, cumulative `r^2`
  subcycling, diffusive reflux, rollback, and serial parity (`0.52.0`)
- [x] transactional owner-only patch-tree `R-T-H-T-R` composition, complete
  bookkeeping synchronization, outer rollback, and serial parity (`0.53.0`)
- [x] rank-local sparse patch payloads, exact owner scatter/gather, and
  same-hierarchy owner-map migration (`0.54.0`)
- [x] direct sparse owner chemistry with distributed average-down, ghost
  refresh, adjacent PPM exchange, rollback, and serial parity (`0.55.0`)
- [x] direct sparse recursive hydro with subcycling, flux registers,
  cross-owner PPM reconciliation, rollback, and serial parity (`0.56.0`)
- [x] direct sparse recursive molecular transport with cumulative `r^2`
  subcycling, diffusive reflux, rollback, and serial parity (`0.57.0`)
- [x] transactional direct sparse `R-T-H-T-R` composition with exact outer
  rollback, call accounting, conservation, and serial parity (`0.58.0`)
- [x] explicit-plan topology-changing sparse regrid with rebuilt ownership,
  overlap retention, rollback, conservation, and serial parity (`0.59.0`)
- [x] tag-driven sparse topology rebuild through four levels with collective
  plan agreement, no-op handling, and rollback (`0.60.0`)
- [x] packed point-to-point same-hierarchy owner migration with one message per
  changed patch and exact payload reconstruction (`0.61.0`)
- [x] point-to-point adjacent sparse sibling halos with exact narrow/PPM
  boundary payloads and cross-owner transfer accounting (`0.62.0`)
- [x] direct sparse child-interior transfer to parent owners for average-down
  and physics synchronization (`0.63.0`)
- [x] direct sparse parent-state fanout once per distinct remote child owner
  for final ghost refresh (`0.64.0`)
- [x] broadcast-free sparse recursive hydro/transport with direct interval,
  boundary-flux, and shared-correction payloads (`0.65.0`)
- [x] replica-free explicit-plan sparse topology rebuild with owner-local
  prolongation and direct old/new overlap transfer (`0.66.0`)
- [x] owner-local solution tagging with compact plan agreement and no
  materialized field tree (`0.67.0`)
- [x] subcycle-weighted deterministic owner assignment and preservation across
  sparse explicit/tag-driven regrids (`0.70.0`)
- [x] owner-local sparse hyperbolic/parabolic timestep evaluation with one
  communicator-wide minimum and no field gather (`0.71.0`)
- [x] namelist-driven sparse MPI AMR time loop, periodic tag regrid, final
  composite diagnostics, and rank-count-invariant CSV output (`0.72.0`)
- [x] versioned sparse patch-tree checkpoint and owner-map-free restart across
  different MPI communicator sizes (`0.73.0`)

AMR parity metrics:

- mass conservation
- energy conservation
- coarse/fine interface behavior
- shock crossing refinement boundaries
- comparison with equivalent-resolution uniform-grid runs

---

## Phase 13 — embedded boundaries

Implement as a separate subsystem after AMR is stable.

Required capabilities:

- [x] geometry representation / nodal level set (`0.74.0`)
- [x] cut-cell detection (`0.74.0`)
- [x] cell-volume fraction (`0.74.0`)
- [x] face-area fraction (`0.74.0`)
- [x] cut-interface length, centroid, and fluid normal (`0.75.0`)
- [x] stationary reactive slip-wall flux and cut-cell source (`0.76.0`)
- [x] Cartesian face-fraction flux divergence (`0.77.0`)
- [x] first-order conservative FluxRedist and state update (`0.78.0`)
- [x] zeroth-order weighted StateRedist with overlapping neighborhoods
  (`0.79.0`)
- [x] piecewise-constant reactive face fluxes and complete EB hydro update
  (`0.80.0`)
- [x] namelist-driven plane/circle EB hydro application, active-cell CFL,
  volume-weighted diagnostics, and geometry-aware CSV (`0.81.0`)
- [x] active-cell chemistry and transactional EB Strang splitting (`0.82.0`)
- [x] active-stencil characteristic PLM and open-face-centroid flux
  interpolation (`0.83.0`)
- [x] fluid-volume centroids and second-order StateRedist neighborhood
  reconstruction, centroid limiting, and recipient bounds (`0.84.0`)
- [x] static aligned two-level EB volume-weighted average-down, composite
  conservation, and reactive EOS transaction (`0.85.0`)
- [x] EB open-face flux register, subcycle accumulation, cut-cell re-reflux,
  fine-recipient transfer, and reactive transaction (`0.86.0`)
- [x] EB PCM prolongation, coarse-time patch-boundary fill, ratio subcycling,
  hydro reflux, and average-down (`0.87.0`)
- [x] input-driven static EB AMR hierarchy, two-level CFL time loop, and
  coarse/fine output (`0.88.0`)
- [x] temperature-gradient tagging, conservative single-patch movement and
  resizing, retained fine overlap, and periodic regrid cadence (`0.89.0`)
- [x] optional fine-patch collapse, root-only advance, and PCM re-creation
  (`0.90.0`)
- [x] lifecycle-aware two-level and root-only EB AMR chemistry composition,
  rollback, and regular-reference parity (`0.91.0`)
- [x] lifecycle-aware serial EB AMR checkpoint, stop, restart, and uninterrupted
  field parity (`0.92.0`)
- [x] deterministic two-level EB multipatch planning, conservative topology
  changes, subcycled hydro/reflux, and Strang chemistry kernel (`0.93.0`)
- [x] input-driven two-level EB multipatch lifecycle, CFL selection, periodic
  regridding, chemistry/hydro advance, and per-child output (`0.94.0`)
- [x] transactional EB multipatch formatted checkpoint/restart with exact
  child topology and public split-run parity (`0.95.0`)
- [x] configured static two-level EB fine patch on an outflow physical
  boundary, including physical-side ghost closure, interface-only reflux, and
  public output qualification (`0.96.0`)
- [x] one-sided physical-boundary tagging and dynamic single-/multipatch
  planning, topology changes, hydro, and output (`0.97.0`)
- [x] static three-level EB composite integration, deepest-to-root
  average-down, reactive EOS recovery, and rollback (`0.98.0`)
- [x] static three-level reactive EB recursive subcycling and reflux with a
  regular finest coarse/fine interface (`0.99.0`)
- [x] EB-cut nested-interface mass, energy, and species conservation closure
  over active middle recipients (`0.100.0`)
- [x] static three-level active-cell chemistry and recursive EB hydro composed
  as a transactional Strang step (`0.101.0`)
- [x] public static three-level hierarchy construction, CFL time loop, and
  per-level output (`0.102.0`)
- [x] dedicated transactional static three-level checkpoint/restart with
  uninterrupted field parity (`0.103.0`)
- [x] tag-driven transactional finest-patch movement and resizing inside a
  fixed middle level (`0.104.0`)
- [x] transactional checkpoint/restart of dynamic three-level finest topology
  and regrid cadence (`0.105.0`)
- [x] single-level EB mixture molecular transport with open-area divergence,
  StateRedist, and symmetric reactive composition (`0.106.0`)
- [x] two-level single-patch EB molecular transport with fine subcycling,
  time-interpolated exterior states, diffusive reflux, and average-down
  (`0.107.0`)
- [x] three-level single-patch EB molecular transport with nested subcycling,
  per-interface diffusive reflux, and deepest-first synchronization (`0.108.0`)
- [x] two-level sibling-patch EB molecular transport with one coarse update,
  per-child subcycling/reflux, set-wide synchronization, and rollback
  (`0.109.0`)
- [x] deterministic MPI ownership of root tiles and sibling EB patches with
  subcycle-weighted work and owner-authoritative synchronization (`0.110.0`)
- [x] owner-only active-cell MPI chemistry for root tiles and sibling EB
  patches with collective commit and rollback (`0.111.0`)
- [x] owner-only MPI EB AMR hydro with one root physics owner, child-owner
  subcycling, flux-register reflux, and rollback (`0.112.0`)
- [x] owner-only MPI EB AMR molecular transport with SSPRK2 stages,
  child-owner diffusive reflux, EB-cut closure, and rollback (`0.113.0`)
- [x] owner-only MPI EB AMR `R-T-H-T-R` composition with outer rollback and
  exact per-operator committed accounting (`0.114.0`)
- [x] sparse rank-local EB root-tile and child payload storage with explicit
  replicated materialization boundary (`0.115.0`)
- [x] direct active-cell chemistry on sparse MPI EB owners with transactional
  materialize/average-down/re-scatter synchronization (`0.116.0`)
- [x] direct sparse child-owner restriction and root-tile-owner reactive
  average-down without complete hierarchy materialization (`0.117.0`)
- [x] sparse-input/output MPI EB AMR full-physics transaction with one central
  replicated `T-H-T` compatibility window (`0.118.0`)
- [x] direct sparse MPI EB AMR hydro with root-level temporary synchronization
  and owner-local child subcycling/reflux (`0.119.0`)
- [x] direct sparse MPI EB AMR SSPRK2 transport with root-level temporary
  synchronization, owner-local child diffusive reflux, and distributed
  cut-interface conservation closure (`0.120.0`)
- [x] end-to-end sparse MPI EB AMR `R-T-H-T-R` composition without a
  replicated fine-child compatibility window (`0.121.0`)
- [x] targeted point-to-point sparse EB child restriction only to intersecting
  root tile owners with exact transfer accounting (`0.122.0`)
- [x] targeted point-to-point direct sparse EB hydro root gather, distinct
  child-owner bundle, correction round trips, and tile scatter (`0.123.0`)
- [x] targeted point-to-point direct sparse EB SSPRK2 transport root traffic,
  final blend scatter, and distributed cut-interface closure (`0.124.0`)
- [x] owner-local sparse EB hydro/transport CFL timestep selection with
  targeted root gather and coarse-interval child scaling (`0.125.0`)
- [x] public sparse EB full-physics time loop with dynamic stable-step
  selection, exact target-time clipping, and committed-prefix accounting
  (`0.126.0`)
- [x] transactional explicit sparse EB topology rebuild with serial overlap
  retention, deterministic owner recomputation, and one-copy post-regrid
  storage (`0.127.0`)
- [x] root-owner temperature tagging and scheduled sparse EB topology rebuild
  inside the public full-physics clock with whole-step rollback (`0.128.0`)
- [x] direct sparse EB regrid restriction, new-owner PCM root assembly, and
  old-owner to new-owner overlap migration with exact traffic accounting
  (`0.129.0`)
- [x] targeted root-only sparse EB materialization for checkpoint/output
  adapters with unallocated non-root fields (`0.130.0`)
- [x] root-selected formatted sparse checkpoint and root/child CSV writers with
  collective I/O status and serial checkpoint compatibility (`0.131.0`)
- [x] root-only formatted checkpoint read and direct root-to-owner sparse
  restart scatter with collective metadata rollback (`0.132.0`)
- [x] geometry-only replicated child topology for direct sparse checkpoint
  restart without replicated child state or temperature (`0.133.0`)
- [x] arbitrary-depth, branching geometry-only EB patch-tree topology with
  transactional whole-tree rebuild (`0.134.0`)
- [x] arbitrary-depth reactive EB numerical hierarchy with conservative
  synchronization and transactional topology/state migration (`0.135.0`)
- [x] finite-halo owner-tiled root hydro for the replicated MPI EB AMR path
  with exact serial parity and bounded work accounting (`0.136.0`)
- [x] point-to-point sparse root halo exchange and owner-tiled EB hydro with
  targeted result routing and exact traffic/work accounting (`0.137.0`)
- [x] zero-gather sparse EB hydro/transport timestep selection directly on
  root tile and child owners (`0.138.0`)
- [x] zero-traffic final sparse SSPRK2 root blend and EOS recovery directly on
  root tile owners (`0.139.0`)
- [x] point-to-point finite-halo sparse SSPRK2 root Euler stages with exact
  target-band accounting and deterministic result routing (`0.140.0`)
- [x] seam-isolated cyclic finite bands for periodic-y sparse root transport
  targets with exact fragment/work accounting (`0.141.0`)
- [x] compact EB flux-register storage and patch-bounded sparse transport
  correction round trips (`0.142.0`)
- locally resolved multilevel EB redistribution and arbitrary-depth physics
  recursion
- [x] first-order isothermal and no-slip embedded-wall transport (`0.180.0`)
- catalytic embedded-wall boundary conditions

Primary PeleC regression references include:

- EB-FlowPastCylinder
- EB-ConvergingNozzle
- EB-TaylorCouette
- EB-BluffBody

---

## Phase 14 — LES

Implement:

- SGS viscosity
- SGS diffusivity
- PeleC-compatible LES models
- dynamic-model support where required

Primary regression case:

- HIT

LES remains separate from the base Navier–Stokes implementation.

---

## Phase 15 — particles and spray

This phase is intentionally late because it depends on stable fluid, transport, MPI, and eventually AMR infrastructure.

Required capabilities:

- Lagrangian particle storage
- particle ownership and MPI migration
- drag
- heat transfer
- evaporation
- breakup
- injection
- two-way coupling
- AMR interaction

Primary PeleC regression references:

- Spray-A-Wbreakup
- Spray-Evaporation
- Spray-Jet
- Spray-Multijet
- Spray-EB

---

## Phase 16 — accelerator support

Only begin after CPU numerical parity is stable.

Order:

1. serial CPU
2. OpenMP CPU
3. MPI
4. MPI + OpenMP
5. accelerator backend

Preferred portable accelerator approaches:

- OpenMP target
- OpenACC

Avoid hard-wiring CUDA Fortran into core physics kernels unless profiling demonstrates a clear need.

---

# I/O strategy

Initial formats:

- CSV for 1D diagnostics
- VTK for simple field output

Later:

- HDF5
- XDMF
- restart/checkpoint files

Checkpoint state must include at least:

- mesh hierarchy
- state arrays
- simulation time
- timestep / iteration count
- AMR metadata

Restart reproducibility should be an automated test.

---

# Parity strategy

Every feature should follow:

```text
PeleC reference case
        ↓
reference data
        ↓
PeleF run
        ↓
automated comparison
```

Possible metrics:

- L1/L2/L∞ field error
- shock position
- contact position
- peak pressure
- peak temperature
- flame position
- ignition delay
- total mass
- total momentum
- total energy
- species conservation

Visual agreement alone is never sufficient for parity acceptance.

---

# Release-level milestones

## PeleF-0

- serial
- single species
- Euler equations
- uniform mesh

## PeleF-1

- PeleC-style Godunov
- multispecies
- NASA thermodynamics

## PeleF-2

- transport
- chemistry
- reacting Navier–Stokes
- PMF-capable uniform-grid solver

## PeleF-3

- MPI
- AMR
- subcycling
- reflux
- dynamic refinement

## PeleF-4

- embedded boundaries
- LES

## PeleF-5

- particles
- spray

## PeleF-6

- GPU / accelerator execution
- large-scale performance tuning

---

# PR strategy

Do not create one enormous translation PR.

Example sequence:

```text
PR 001  project skeleton + architecture
PR 002  state representation
PR 003  ideal-gas EOS
PR 004  primitive/conserved conversion
PR 005  Rusanov solver
PR 006  finite-volume update
PR 007  Sod regression test
PR 008  CFL timestep
PR 009  PLM
PR 010  PeleC-style Riemann solver
...
```

## 0.143.0 sparse MPI compact child transport-context gates

Replace the complete root transport bundle sent to every distinct child owner
with a per-child four-edge start/end exterior context and the compact coarse
flux-register mismatch. Advance fine transport only on the child owner, return
the evolved fine field and accumulated register to the root physics owner,
apply reflux there in deterministic child order, and return only the corrected
fine field.

Require remote children to allocate no complete root state, temperature, or
x/y-flux field. Verify exact complete-root/context exterior parity, a strict
context payload reduction, three messages per remote child per Euler stage,
and unchanged numerical, work, rollback, and public-clock gates at one, two,
four, and eight ranks before the complete serial regression.

## 0.144.0 sparse MPI compact child-local reflux gates

Generalize EB reactive reflux to operate on a globally indexed coarse support
rectangle containing the patch-plus-two footprint. Send that support with the
fine exterior/register context, execute fine subcycling and reactive reflux on
the child owner, keep the corrected fine field there, and return only corrected
coarse support for deterministic root-owner merge.

Require support/full reflux parity, a combined context payload smaller than
the former root bundle, exactly two messages per remote child per Euler stage,
and unchanged numerical, owner-work, rollback, and public-clock gates at one,
two, four, and eight ranks before the complete serial regression.

## 0.145.0 compact coarse interface-flux gates

Generalize coarse EB flux-register accumulation to accept globally indexed
x/y face rectangles that contain the active coarse/fine interface. Preserve
the complete-root entrypoint as a wrapper and switch sparse MPI child-register
initialization to the compact call.

Require compact/full bitwise correction parity, strict payload reduction,
transactional rejection of incomplete support, and unchanged sparse numerical,
work, traffic, clock, scheduled-regrid, and rollback gates at one, two, four,
and eight ranks before the complete serial regression.

## 0.146.0 direct root-tile coarse-flux routing gates

Retain globally indexed x-flux rows and uniquely owned y-faces on every root
tile owner. Route only child-intersecting fragments to each child owner,
assemble the compact interface rectangles there, and initialize the coarse
flux register without root-physics-owner involvement. Remove the register from
the root-to-child state context.

Require complete fragment coverage, combined compact payload reduction, exact
distribution-derived message counts, and unchanged sparse numerical, owner-
work, limiter, clock, scheduled-regrid, and rollback gates at one, two, four,
and eight ranks before the complete serial regression.

## 0.147.0 compact exterior state-context gates

Generalize reactive child exterior-context extraction to globally indexed
coarse start/end state and temperature support containing the patch-plus-one
footprint. Preserve complete-root extraction as a wrapper.

Require strictly smaller support, bitwise complete/support exterior parity,
transactional invalid-support rejection, and unchanged serial and sparse MPI
reactive EB AMR regressions before direct root-tile state routing.

## 0.148.0 direct root-tile state/support routing gates

Retain start, uncorrected-end, and current corrected state/temperature on each
root transport tile. Route patch-plus-two fragments directly to child owners,
extract exterior context there, and return child-local reflux corrections to
the intersecting tile owners in deterministic child order.

Remove the root-owner child context/correction path and final corrected-row
scatter. Require reduced compact payload, exact tile/child traffic counts,
serial numerical parity, owner work, limiter, clock, scheduled-regrid,
cut-boundary conservation, and rollback gates at one, two, four, and eight
ranks before the complete serial regression.

## 0.149.0 owner-local root transport result gates

Remove the complete owned transport result sent from every remote root tile to
the root physics owner after each Euler stage. Retain stage state and fluxes
only in owner-local tile records used by direct child routing and final commit.

Replace complete-root cut-boundary flux inspection with tile-local physical-
boundary contributions and one communicator-wide `nvar` sum. Require exact
halo/child-only point-to-point traffic, root-only cyclic parity, cut-interface
conservation, owner work, limiter, clock, scheduled-regrid, and rollback gates
at one, two, four, and eight ranks before the complete serial regression.

## 0.150.0 compact sparse hydro child-context gates

Replace the complete hydro root bundle sent to every distinct remote child
owner with one child-specific message containing four-edge start/end context,
current patch-plus-two corrected state/temperature, and intersecting coarse
x/y flux support. Run coarse-register accumulation and reactive reflux on the
child through existing globally indexed support APIs, then return only the
corrected support.

Require the compact payload to be strictly smaller than the former root bundle
and remove the distinct-owner bundle allowance from exact traffic. Preserve
serial field parity, deterministic child order, owner work, conservation,
clock, scheduled-regrid, and rollback gates at one, two, four, and eight ranks
before the complete serial regression.

## 0.151.0 direct hydro coarse-flux routing gates

Keep hydro x rows and uniquely owned y-faces on root tile owners. Send only
state and temperature in remote tile results and in the compact root-to-child
context, then route every intersecting x/y fragment directly to the child
owner for covered, finite interface-register assembly.

Require the root physics owner to allocate no complete hydro flux field and
derive exact traffic from remote halo, state-result, final-scatter, child
state/correction, and tile/child flux-intersection messages. Preserve serial
field parity, deterministic child order, owner work, conservation, clock,
scheduled-regrid, and rollback gates at one, two, four, and eight ranks before
the complete serial regression.

## 0.152.0 owner-local hydro result gates

Retain start, uncorrected-end, and current corrected hydro state/temperature on
root tile owners. Route patch-plus-two state directly to child owners, extract
the four-edge context there, and return reflux corrections directly to every
intersecting tile owner in deterministic child order.

Remove remote tile results, complete root hydro result allocation, root-owner
support merge, and final row scatter. Require exact traffic to contain only
finite-band halos and direct state, flux, and correction fragments. Preserve
serial field parity, owner work, conservation, clock, scheduled-regrid, and
rollback gates at one, two, four, and eight ranks before the complete serial
regression.

## 0.153.0 arbitrary-depth reactive EB patch-tree timestep gates

Move the qualified active-cell EB CFL calculation into a shared numerical
module while retaining the existing driver API. Traverse every runtime tree
node and reduce its local limit into root time with the cumulative product of
all ancestor refinement ratios. Do not add a fixed level bound or materialize a
second hierarchy.

Require independent four-level branching parity, a deepest-node limiting case,
finite positive output, deterministic zero on invalid input, and a read-only
state/temperature contract. Run the complete serial Debug and Release suite and
the existing OpenMPI one-, two-, four-, and eight-rank gates before acceptance.

## 0.154.0 arbitrary-depth reactive EB patch-tree hydro gates

Generalize the qualified fixed-depth EB advance into a node-recursive
transaction. Retain parent start/end state for child time interpolation,
subcycle each relation by its runtime ratio, and own one flux register per
ordered child. Reflux, average down, and close each refined subtree against its
outer-boundary flux before a deepest-first final synchronization and atomic
commit.

Require a four-level branching schedule and composite density, energy, and
species conservation; compare a runtime three-level chain with the qualified
fixed-depth implementation; and verify exact rollback plus zero counts after a
late recursive rejection. Run the complete serial Debug and Release suite and
the existing OpenMPI one-, two-, four-, and eight-rank gates before acceptance.

## 0.155.0 arbitrary-depth reactive EB patch-tree chemistry gates

Traverse every runtime patch with its own EB active mask and the qualified 2D
chemistry integrator. Expose standalone chemistry plus one transactional
`chemistry(dt/2) -> recursive hydro(dt) -> chemistry(dt/2)` operation. Keep the
accepted tree and public chemistry/hydro counts unchanged until final
deepest-first synchronization and validation succeed.

Require exact four-level branching chemistry and hydro schedules, composite
mass/energy/species closure, reaction activity, and rollback after a valid
chemistry prefix followed by rejected hydrodynamics. Compare a runtime
three-level chain with the established fixed-depth Strang implementation. Run
the complete serial Debug and Release suite and the existing OpenMPI one-,
two-, four-, and eight-rank gates before acceptance.

## 0.156.0 arbitrary-depth reactive EB patch-tree transport gates

Generalize the qualified fixed-depth EB transport Euler stage into a
node-recursive operation with parent start/end interpolation, ratio subcycling,
one diffusive flux register per child, subtree conservation closure, and
deepest-first synchronization. Compose two complete Euler trees with node-wise
SSPRK2 blending and EOS temperature recovery on one private candidate.

Require fixed three-level field/temperature parity and exact recursive counts,
then exercise a four-level branching topology with composite conservation,
positive limiter theta, changed state, thermodynamic validity, and exact
rollback. Run the complete serial Debug and Release suite and the existing
OpenMPI one-, two-, four-, and eight-rank gates before acceptance.

## 0.157.0 arbitrary-depth reactive EB patch-tree full-physics gates

Compose the qualified all-node reaction traversal, recursive SSPRK2 transport,
and recursive hydrodynamics as one private `R-T-H-T-R` candidate. Defer state,
temperature, minimum transport theta, and all three per-level count vectors
until final deepest-first synchronization and tree validation succeed.

Require fixed three-level full-physics field/temperature parity and exact
chemistry, transport, and hydro schedules. Exercise the four-level branching
tree for actual recursive counts, composite conservation, positive finite
thermodynamics, and rollback after valid reaction and transport prefixes. Run
the complete serial Debug and Release suite and the existing OpenMPI one-,
two-, four-, and eight-rank gates before acceptance.

## 0.158.0 public arbitrary-depth reactive EB patch-tree time-loop gates

Add an all-node stable-step selector that reduces both hyperbolic and active
explicit transport limits after cumulative refinement scaling. Compose it
with the qualified `R-T-H-T-R` transaction in a caller-owned target-time loop
with exact stop clipping and a maximum-step bound.

Treat each interval as one candidate transaction and publish tree state,
clock, step count, minimum accepted interval, limiter minimum, and accumulated
per-level physics counts only after acceptance. Require independent two-step
parity, exact accounting, committed-prefix behavior at the step bound, and
exact first-step rollback. Run the complete serial Debug and Release suite and
the existing OpenMPI one-, two-, four-, and eight-rank gates before acceptance.

## 0.159.0 MPI arbitrary-depth EB patch-tree ownership gates

Introduce a topology-matched distribution with one owner per runtime node,
deterministic greedy placement, cumulative-subcycle work weighting, and exact
per-rank accounting. Retain a replicated numerical tree initially and publish
owner-authoritative state and temperature through one private all-rank
candidate.

Require collective topology/control consensus, exact four-level branching
owner accounting, all-rank field parity, rank-local invalid-state rollback,
and inconsistent-control rejection at one, two, four, and eight ranks. Run the
complete existing MPI and serial Debug/Release suite before acceptance. Defer
sparse nonowner storage, direct migration, and owner-local physics routing.

## 0.160.0 MPI sparse arbitrary-depth EB patch-tree storage gates

Keep topology and ownership replicated while allocating each node's numerical
fields only on its owner. Add an explicit replicated materialization boundary
for compatibility and direct point-to-point migration from each old owner to
the corresponding new owner through one private sparse candidate.

Require exact owner-local allocation accounting, pre/post-migration field
parity, exact changed-owner transfer accounting, and collective invalid-map
rollback at one, two, four, and eight ranks. Run the complete existing MPI and
serial Debug/Release suite before acceptance. Defer distributed sparse
timestep reduction and owner-local recursive physics routing.

## 0.161.0 MPI owner-local arbitrary-depth EB timestep gates

Evaluate hyperbolic and explicit-transport stability limits only on each
node's sparse owner. Convert every local interval to root time with the
cumulative refinement product, reduce the communicator minimum, and publish
exact active-node accounting without materializing a complete tree.

Require exact serial-selector parity after owner rotation and collective
rejection of rank-local CFL disagreement at one, two, four, and eight ranks.
Run the complete existing MPI and serial Debug/Release suite before acceptance.
Defer owner-local recursive hydro, transport, chemistry, and clock routing.

## 0.162.0 MPI owner-local arbitrary-depth EB chemistry gates

Advance chemistry only on each sparse node owner. Recover local temperatures,
then synchronize relations deepest-first. Copy shared-owner children locally
and send conserved state once from a distinct child owner to the parent owner
before applying the serial EB average-down kernel.

Require exact serial state/temperature parity, per-level chemistry counts,
map-derived direct-transfer counts, and collective control-mismatch rollback at
one, two, four, and eight ranks. Run the complete existing MPI and serial
Debug/Release suite before acceptance. Defer owner-local recursive hydro,
transport, and public clock routing.

## 0.163.0 MPI sparse arbitrary-depth EB composite-integral gates

Expose whole-tree and selected-subtree conserved integrals directly on sparse
owners. Traverse the replicated topology recursively, mask direct-child
coverage before integrating a parent, and reduce only owner-local contributions
without materializing fields.

Require serial parity for the complete four-level branching tree and every
subtree, exact topology-derived contributing-node counts, and collective
selector-disagreement rejection with neutral outputs at one, two, four, and
eight ranks. Run the complete existing MPI and serial Debug/Release suite
before acceptance. Use this reduction as the conservation prerequisite for
owner-local recursive hydro; defer transport and public clock routing.

## 0.164.0 MPI owner-local arbitrary-depth EB hydro gates

Route the established recursive hydro schedule directly over sparse owners.
Send compact parent-time exterior context to a distinct child owner, return
fine fluxes to the parent register after every substep, and route corrected
child state through parent-owner reflux and ordered average-down. Reuse sparse
subtree reductions for conservation closure and commit one private candidate.

Require serial field/composite parity, exact per-level advance scheduling,
topology/map-derived grouped-transfer counts, and collective interval-mismatch
rollback at one, two, four, and eight ranks. Run the complete existing MPI and
serial Debug/Release suite before acceptance. Defer owner-local recursive
transport and the public sparse full-physics clock.

## 0.165.0 MPI owner-local arbitrary-depth EB transport gates

Route both SSPRK2 Euler stages over sparse node owners. Reuse compact parent-
time context, direct fine-flux return, parent-owner registers, child-state
reflux round trips, ordered average-down, and owner-local subtree conservation
closure. Blend the accepted and second-stage fields and recover EOS temperature
only on the node owner, then perform one direct deepest-first synchronization.

Require serial field, temperature, limiter, and composite-integral parity;
exact per-level Euler scheduling; topology/map-derived grouped-transfer counts;
and collective interval-mismatch rollback at one, two, four, and eight ranks.
Run the complete MPI and serial Debug/Release suite before acceptance. Defer
owner-local full-physics composition and the public sparse clock.

## 0.166.0 MPI owner-local arbitrary-depth EB full-physics gates

Compose the qualified sparse chemistry, SSPRK2 transport, and recursive hydro
entrypoints as one private `R-T-H-T-R` candidate. Establish outer consensus
before optional-physics branching, accumulate each operator and transfer class
separately, reduce the two transport limiter minima, and publish only after
final sparse validation.

Require exact serial split scheduling and topology/map-derived traffic, plus
field, temperature, limiter, and composite-integral parity. Reject a rank-local
interval mismatch before mutation with zero public diagnostics at one, two,
four, and eight ranks. Run the complete MPI and serial Debug/Release suite
before acceptance. Defer the owner-local target-time clock.

## 0.167.0 MPI owner-local arbitrary-depth EB clock gates

Wrap the owner-local timestep selector and sparse `R-T-H-T-R` transaction in a
public target-time loop. Establish clock/control consensus before the loop,
recompute a stable interval before every attempted step, clip the final interval
to the requested target, and publish time, step, minima, advances, and transfer
counts only with each accepted sparse candidate. Preserve the committed prefix
on a later failure or step ceiling.

Require exact one-step target-time and minimum-dt parity, serial field/limiter/
composite parity, exact timestep-node and physics/traffic accounting, clock-
control mismatch rollback, and maximum-step rollback at one, two, four, and
eight ranks. Run the complete MPI and serial Debug/Release suite before
acceptance. Defer arbitrary-depth dynamic tagging and checkpoint I/O.

## 0.168.0 serial arbitrary-depth EB tagged-rebuild gates

Synchronize the accepted numerical tree deepest-first, plan temperature tags
for every prospective parent independently, and build deterministic
parent-major child plans through the configured maximum depth. Keep EB geometry
construction behind a caller callback, and treat parents below the tagger's
minimum stencil extent as terminal branches.

Require a root-only hot-cell case to reach three levels, preserve the complete
composite conserved vector, and retain exact parent ownership. Repeating the
same plan must be a field-exact no-op. A rejecting geometry callback must leave
the complete accepted tree unchanged, and removing every temperature gradient
must collapse the tree to its root while preserving the composite integral.
Run the complete serial Debug and Release suite before acceptance. Defer
owner-local MPI tag planning/migration and arbitrary-depth checkpoint I/O.

## 0.169.0 MPI owner-local arbitrary-depth EB tagged-rebuild gates

Evaluate each prospective parent only on its owner, reduce compact tag-plan
metadata, rebuild caller-defined EB geometry collectively, and assign the
candidate topology with the deterministic work model. Initialize and migrate
candidate fields through direct parent/child and old-owner/new-owner traffic;
never materialize a complete numerical tree.

Require root-only creation through three levels, changed-plan retained overlap,
an unchanged-plan exact no-op, and a tag-free collapse. Compare plan metadata,
fields, temperatures, and composite integrals with the serial reference at one,
two, four, and eight ranks. Require exact topology-derived transfer accounting
and collective rollback on rank-dependent or invalid criteria. Run the complete
MPI and serial Debug/Release suite before acceptance. Defer arbitrary-depth
checkpoint/restart and composite output.

## 0.170.0 serial arbitrary-depth EB checkpoint gates

Store the complete ordered branching topology, every root/child EB metric,
conserved state, temperature, species order, and lifecycle counter in a
distinct versioned formatted stream. Reconstruct a private candidate and
recover general-EOS temperature before publication.

Require a four-level branching round trip with topology, field, temperature,
and metadata parity. Reject an insufficient configured depth, exchanged species
order, invalid lifecycle metadata, malformed topology, invalid EOS state, or a
missing terminal marker without publishing partial state. Run the complete
serial Debug/Release suite before acceptance. Defer sparse MPI root-only I/O,
rank-neutral restart, and composite output.

## 0.171.0 sparse MPI arbitrary-depth EB checkpoint gates

Gather each numerical node only to a selected I/O root, write the qualified
serial tree format there, and store no owner map. On restart, read only on the
root, broadcast compact topology/geometry, recompute the current deterministic
distribution, and scatter fields directly to new owners.

Require exact topology-derived gather/scatter traffic and serial field,
temperature, and lifecycle metadata parity at one, two, four, and eight ranks.
Change the work exponent across restart to prove redistribution. Reject
rank-dependent controls and species identity before traffic with neutral public
outputs. Run the complete MPI and serial Debug/Release suite before acceptance.
Defer arbitrary-depth composite output.

## 0.172.0 arbitrary-depth EB composite output gates

Write one serial CSV from the complete branching tree. Exclude cells covered
by direct children and include level/patch identity, EB geometry diagnostics,
conserved and primitive fields, temperature, and ordered mass fractions.

For sparse MPI, gather each node directly to a caller-selected writer root and
invoke the same serial writer only there. Require exact remote-node transfer
counts and exact topology-derived leaf-row counts at one, two, four, and eight
ranks. Run the complete MPI and serial Debug/Release suite before acceptance.
Defer the runnable arbitrary-depth 2D EB application lifecycle.

## 0.173.0 runnable serial arbitrary-depth EB application gates

Add a dedicated executable that reuses the established reactive 2D, EB, and
AMR namelists while selecting the arbitrary-depth numerical tree. Initialize or
restart the tree, apply scheduled recursive tags, select stable all-node root
steps, advance `R-T-H-T-R`, invoke scheduled/final checkpoints, and write one
composite CSV.

Require a public input case to populate four levels and finish with valid EB,
thermodynamic, species-closure, spacing, identity, and time data. Run the full
210-test serial suite in GNU Fortran Release and bounds/FPE-checked Debug.
Defer application-level checkpoint/restart comparison and public sparse MPI
application integration.

## 0.174.0 public patch-tree application restart gates

Exercise the installed serial arbitrary-depth EB executable through three
separate processes: uninterrupted reference, checkpoint-stop after the first
committed step, and continuation from that checkpoint. Keep recursive initial
tagging and periodic regridding enabled so the file carries a four-level
numerical hierarchy and the resumed clock retains its global cadence.

Require a structurally complete checkpoint, an intermediate stopped output,
exact final times, identical final composite identities and columns, and
bounded numeric differences for every field. Run the complete 214-test suite
in GNU Fortran Release and bounds/FPE-checked Debug. Defer public sparse-MPI
application integration and explicit input/checkpoint compatibility hashes.

## 0.175.0 public sparse-MPI application gates

Connect the existing sparse arbitrary-depth EB APIs to an installed executable
that reuses the serial application inputs. Cover fresh initialization,
recursive regridding, full physics, checkpoint/restart calls, collective
integrals, and selected-root output without materializing child fields on
non-owners. Compare 1/2/4/8-rank composite fields and retain all serial and MPI
Debug/Release gates. Defer cross-rank application restart composition and
elimination of the temporary replicated root startup field.

## 0.176.0 public sparse-MPI cross-rank restart gates

Exercise the installed sparse-MPI arbitrary-depth EB executable through four
separate processes: an uninterrupted one-rank reference, a two-rank
checkpoint-stop run, and independent four- and eight-rank continuations from
that file. Change the MPI work exponent across the checkpoint boundary so the
test covers ownership-policy recomputation in addition to communicator-size
redistribution.

Require a structurally complete four-level checkpoint, an intermediate stop,
exact final times, identical final composite identities and columns, and
bounded differences for every numeric field. Run the complete MPI gate chain
and all 214 serial tests in GNU Fortran Release and bounds/FPE-checked Debug.
Defer elimination of the temporary replicated root startup field.

## 0.177.0 owner-local public sparse-MPI startup gates

Construct the fresh root topology and distribution before numerical state.
Run the established reactive initializer only on the root-node owner and add a
collective root-only sparse initializer that rejects any non-owner numerical
allocation. Move the owner arrays into the sparse node and require the source
arrays to be unallocated before initial recursive regridding.

Exercise the resulting startup path at one, two, four, and eight ranks and
retain its complete-field parity, cross-rank checkpoint/restart parity, the
complete MPI gate chain, and all 214 serial tests in GNU Fortran Release and
bounds/FPE-checked Debug. Defer explicit application/checkpoint compatibility
fingerprints.

## 0.178.0 public patch-tree checkpoint fingerprint gates

Add a structured compatibility fingerprint to public serial and sparse-MPI
arbitrary-depth EB checkpoints. Include mesh/domain, EB construction, physics,
StateRedist, hierarchy, and regridding controls; exclude continuation length,
output/checkpoint scheduling, communicator size, and ownership weighting.

Require valid serial and 2-to-4/eight-rank MPI restarts to retain complete-field
parity. Require a CFL-mismatched restart to fail transactionally in serial and
MPI, while retaining low-level schema-1 compatibility and all 215 serial and
MPI Debug/Release gates.

## 0.179.0 interface-local multilevel EB closure gates

Replace parent-wide residual spreading after EB reflux and average-down with a
topology-derived local support. For every direct child, mark active unrefined
parent cells in clipped three-by-three neighborhoods of the coarse/fine
interface and normalize the conserved correction by their fluid volume.

Require fixed-depth, multipatch, arbitrary-depth, serial, and sparse-MPI paths
to use the same support rule. Retain density, total-energy, species, EOS,
rollback, and complete-field parity gates in all 215 serial tests and the full
one-, two-, four-, and eight-rank MPI Debug/Release chain. Do not claim exact
AMReX per-neighborhood transfer parity.

## 0.180.0 embedded-wall thermal and viscous gates

Add one validated embedded-wall record to the shared 2D boundary set and keep
its default adiabatic, free slip, and species impermeable. Evaluate first-order
normal Fourier heat transfer for isothermal walls and Newtonian traction plus
moving-wall work for no-slip walls. Insert only the wall-length-weighted flux
into EB cut-cell transport divergence so every existing AMR and sparse-MPI
transport route reuses it without a new stepping interface.

Require direct heat/traction/work signs, exact slip and impermeability zeros,
invalid-distance rejection, cut-cell-only right-hand-side changes, cross-rank
boundary consensus, all 215 serial tests, and the full MPI Debug/Release chain.
Defer catalytic species fluxes, higher-order wall stencils, namelist exposure,
and nondefault-wall checkpoint fingerprints.

## 0.181.0 public single-level embedded-wall controls

Expose wall kind, thermal mode, temperature, and velocity in the established
`&embedded_boundary` namelist. Validate their transport dependencies and apply
them transactionally to the shared boundary set before the single-level public
clock starts. Exercise the installed application with a hot, tangentially
moving no-slip wall and require cut-cell heat and momentum response.

Reject active nondefault values in checkpoint-capable AMR application
preflight until they join every restart compatibility record. Retain all 215
serial tests and the full MPI Debug/Release chain. Defer AMR namelist/checkpoint
exposure, catalytic fluxes, and higher-order wall stencils.

## 0.182.0 restart-safe AMR embedded-wall controls

Reuse the configured boundary builder in fixed-depth, arbitrary-depth, and
sparse-MPI public AMR paths. Advance the fixed-depth checkpoint schemas and the
arbitrary-depth serial/MPI fingerprint with every wall and molecular-transport
compatibility control. Qualify active isothermal-wall checkpoint/restart and
mismatch rollback while preserving the established low-level stepping APIs.

## 0.183.0 EB-safe limited-linear AMR prolongation

Add a conservative MC-limited coarse-to-fine initializer for reactive 2D EB
patches. Apply nonzero slopes only when the parent and all of its children are
regular, retain PCM across EB topology changes, recover child temperatures
through the EOS, and retry a parent with PCM if a linear child is inadmissible.
Qualify analytic linear reproduction, prolong/restrict conservation, cut-parent
fallback, and neutral-output rejection before exposing runtime selection.

Each implementation PR should normally contain:

1. implementation
2. unit tests
3. parity/regression test
4. documentation update

---

# Anti-patterns

Do not use the following workflow:

```text
translate all C++ files automatically
↓
fix compiler errors
↓
assume correctness
```

Instead use:

```text
identify one PeleC capability
↓
understand numerical behavior
↓
implement in Fortran
↓
unit test
↓
compare against PeleC/reference solution
↓
accept parity
↓
move to next capability
```

---

# First implementation batch

The first code batch should contain:

```text
docs/architecture.md
docs/pelec_mapping.md
docs/parity_strategy.md
CMakeLists.txt
src/core/precision_mod.F90
src/core/state_indices_mod.F90
src/physics/eos_ideal_mod.F90
src/hydro/state_conversion_mod.F90
src/hydro/riemann_rusanov_mod.F90
src/hydro/finite_volume_mod.F90
src/driver/time_integrator_mod.F90
app/pelef.F90
tests/unit/
cases/sod/
tools/compare_sod.py
.github/workflows/
```

First hard acceptance criterion:

> Build successfully with gfortran, run a standalone Sod shock tube simulation to completion, preserve positive density/pressure, satisfy conservation checks, and pass automated reference comparison in CI.

Only after that gate passes should PeleC-specific Godunov reconstruction and larger capabilities be added.
