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

Milestone `0.221.0` implements the first bounded form of this strategy for the
configure-time selected serial 0D application: pinned SUNDIALS 7.2.0 CVODE
BDF, serial N_Vector, dense linear algebra, reduced mass fractions, fixed-energy
temperature recovery, and a semi-analytic Jacobian. The dependency is optional
and isolated, requires the official static Fortran-module targets, and rejects
shared-only prefixes to keep the installed executable self-contained with
respect to SUNDIALS. The native adaptive backward-Euler solver remains the
default.
Milestone `0.222.0` adds an opaque handle and generation-checked registry for
up to 64 live contexts. Their calls may be interleaved sequentially with
independent failure and cumulative-work state. Threaded/reentrant execution,
sparse linear algebra, CFD reaction-source dispatch, detailed mechanisms, and
performance qualification remain Phase-7 work.

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

## 0.184.0 fixed-depth public prolongation selection

Expose `pcm` and `linear` in `&eb_amr`, keep PCM as the compatibility default,
and route the selection through static and dynamic two-level, sibling-patch,
and three-level initialization. Qualify the installed hot-wall transport case
with linear prolongation and retain explicit rollback for an unknown method.

Until the selection participates in every restart identity, reject linear
fixed-depth checkpoint/restart requests and keep arbitrary-depth serial and
sparse-MPI patch trees PCM-only. Track the method in those schemas before
lifting either restriction.

## 0.185.0 restart-safe arbitrary-depth prolongation selection

Advance every fixed-depth checkpoint schema and the shared serial/sparse
patch-tree fingerprint with `prolongation_method`. Remove the temporary
PCM-only lifecycle guards only after reads reject a method mismatch before
publication.

Propagate the selected dispatcher through iterative serial tag planning and
rebuild, plus sparse-MPI owner-local candidate construction and final rebuild.
Require collective method agreement and retain direct parent-owner to
child-owner routing. Qualify public fixed-depth, four-level serial, 1/2/4/8
rank, and cross-rank restart cases with linear prolongation.

## 0.186.0 conservative limited-linear cut-parent prolongation

Replace the cut-parent PCM-only branch with slopes measured between active
coarse fluid centroids. Use MC slopes with two-sided support and one-sided
derivatives at the embedded boundary. Reconstruct at fine fluid centroids,
remove their volume-fraction-weighted mean offset, and limit every conserved
component to the active coarse-neighbor envelope.

Retain PCM for covered or other topology-mismatched parents and retain the
existing parent-local PCM retry after failed EOS recovery. Qualify nonconstant
cut-parent children, component bounds, regular-parent exactness, and EB
prolongation/average-down conservation before using the shared dispatcher in
all existing serial and sparse-MPI lifecycles.

## 0.187.0 multidimensional cut-parent prolongation

Replace independent coordinate secants at EB-cut parents with a connected
3-by-3 least-squares fit over coarse fluid-volume centroids. Admit a diagonal
sample only when an open two-face path reaches it. Solve a full-rank
two-dimensional normal matrix, retain the minimum-norm tangent for rank-one
support, and use zero slope when no direction is resolved.

Limit connected coarse-centroid predictions and active fine children to the
local active component envelope. Retain the fine-volume-weighted zero-mean
offset, average-down conservation, EOS recovery, and parent-local PCM retry.
Qualify an interface-tangential affine field defined at fluid centroids before
retaining every fixed-depth, arbitrary-depth, serial, and sparse-MPI lifecycle
gate.

## 0.188.0 rank-recovering cut-parent stencil

Keep the connected 3-by-3 least-squares fit as the compact path. Before
accepting a rank-one or empty result, rebuild the normal system over a 5-by-5
box using a bounded open-face flood fill from the cut parent. Solve the grown
system in two dimensions when it becomes full rank and use its connected
component envelope for limiting.

Qualify a turning fluid path whose compact neighborhood resolves only one
direction but whose grown neighborhood resolves both components of an affine
field. Retain exact fine-volume-weighted average-down, EOS recovery, PCM retry,
and all shared fixed-depth, arbitrary-depth, serial, and sparse-MPI gates.

## 0.189.0 transactional fixed three-level parent regrid

Plan a replacement root-to-middle rectangle from root temperature tags. Before
moving it, restrict the accepted finest contribution into the old middle, then
reuse the two-level regrid transaction to conserve into the root and retain
same-resolution middle overlap. Plan and prolong a replacement finest patch
inside the rebuilt middle's two-cell safety margin.

Publish all three state/temperature fields, both refined geometries, and both
patch descriptors only after EOS validation and a before/after three-level
composite conservation check. Qualify changed parent bounds, valid nested
topology, conservation, and invalid-control rollback. Keep public scheduling
and fixed-depth checkpoint topology for the following lifecycle milestone.

## 0.190.0 public fixed three-level parent lifecycle

Add an explicit fixed-depth parent-regridding control without changing legacy
dynamic-finest inputs. Route both initialization and periodic regrid events
through a parent-first dispatcher and count one committed hierarchy change.
Reject configurations whose minimum root-tagged patch cannot support the
finest two-cell margin.

Advance the dynamic three-level checkpoint schema, persist the policy and both
actual patch descriptors, and reconstruct the stored parent before validating
the finest descriptor. Qualify a public moved-parent hotspot plus uninterrupted
versus checkpoint/restart field parity with topology-derived output sizes.

## 0.191.0 arbitrary-depth outflow-boundary children

Qualify the shared domain-inclusive temperature planner in the public serial
and sparse-MPI arbitrary-depth applications. Drive a hotspot against one
outflow side through four populated levels, copy the current fine boundary
state on that physical side, and retain coarse-time interpolation elsewhere.

Require physical-side flux-register omission, complete composite EB classes,
positive reactive fields, species closure, and exact 1/2/4/8-rank composite
output parity. Leave non-outflow physical boundaries and periodic-seam
children explicit for later work.

## 0.192.0 boundary-touching patch-tree restart

Move the public arbitrary-depth checkpoint/restart hotspot to the qualified
x-upper boundary. Require the uninterrupted, stopped, and restarted composite
outputs to retain four levels and exact physical-side contact.

Reuse checkpoint schema 4 because its stored child bounds and geometry already
identify domain-inclusive patches. Qualify an independent serial restart and a
two-rank checkpoint restarted at four and eight ranks with changed ownership
weight, while retaining fingerprint mismatch rejection and full field parity.

## 0.193.0 boundary-touching arbitrary-depth transport

Enable thermal conduction in both public x-upper patch-tree cases. Exercise
the recursive `r^2` transport schedule, coarse-time child context, physical-
side diffusive register omission, reflux, and composite synchronization in
serial and sparse-MPI execution.

Carry the active transport fingerprint through checkpoint-stop and independent
restart. Retain exact 1/2/4/8-rank fresh-run parity and the two-rank checkpoint
to four-/eight-rank continuation gate without adding duplicate regressions.

## 0.194.0 boundary-touching mixture transport

Enable viscosity, mixture-averaged species diffusion, and barodiffusion beside
Fourier conduction in both public x-upper trees. Exercise correction velocity,
species enthalpy flux, recursive diffusive registers, physical-side omission,
reflux, and composite synchronization as one transport transaction.

Reuse the existing fresh and restart regressions so the complete transport
combination must preserve serial and 1/2/4/8-rank field identity, including a
two-rank checkpoint continued at four and eight ranks under changed ownership.

## 0.195.0 reacting boundary-touching full physics

Enable elementary chemistry in both public x-upper cases while retaining the
complete transport combination. Qualify the full transactional `R-T-H-T-R`
sequence on four populated levels with periodic recursive regridding.

Reuse the fresh and split-run gates so reaction-modified states must retain
serial and 1/2/4/8-rank field identity, plus two-rank checkpoint continuation
at four and eight ranks under a changed ownership weight.

## 0.196.0 restart-persistent transport limiter history

Advance the arbitrary-depth checkpoint envelope to base schema 2 and public
fingerprinted schema 5. Store the cumulative minimum transport limiter theta
with time, step, regrid, and minimum-timestep metadata.

Restore that value before serial or sparse-MPI continuation, including the
selected-root broadcast and changed-rank ownership rebuild. Extend serial and
MPI round-trip gates with a nonneutral synthetic value and keep failed-read
rollback at neutral theta `1`.

## 0.197.0 restart-persistent conservation baseline

Advance the arbitrary-depth checkpoint envelope to base schema 3 and public
fingerprinted schema 6. Store the original `nvar`-component composite integral
after the clock metadata, validate finite values, and restore it before any
continuation output or final diagnostic.

Require exact direct serial and selected-root sparse round trips, collective
rejection of rank-disagreed baselines, and an unallocated optional output on
failed reads. Extend the public two-to-four/eight-rank restart checker to
compare the cumulative conservation error and limiter minimum with the
uninterrupted reference.

## 0.198.0 restart-persistent AMR operator counters

Advance the arbitrary-depth checkpoint envelope to base schema 4 and public
fingerprinted schema 7. Store fixed-capacity, per-level cumulative chemistry,
transport, and hydro patch-advance vectors. Size public-driver vectors to the
configured maximum tree depth so a temporary topology shrink does not discard
deeper-level history.

Require direct serial recovery with capacity beyond the populated depth,
sparse rank consensus, selected-root broadcast, and unallocated optional
outputs after failed reads. Capture the public uninterrupted and restarted
logs and require exact counter-vector equality across the serial and
two-to-four/eight-rank continuation chains.

## 0.199.0 restart-persistent AMR regrid history

Advance the arbitrary-depth checkpoint envelope to base schema 5 and public
fingerprinted schema 8. Store cumulative successful tag/regrid evaluations and
the sum of globally tagged cells. Count an evaluation only after its complete
transaction succeeds, including evaluations that retain the existing
topology.

Require exact direct serial recovery, sparse rank consensus, selected-root
broadcast, and neutral outputs after rejected reads. Extend serial and
two-to-four/eight-rank restart log comparison so both adaptation diagnostics
match the uninterrupted logical run exactly.

## 0.200.0 public branching patch-tree lifecycle

Replace the single-feature public boundary and restart cases with separated
boundary-touching and interior temperature features. Require at least two
leaf-visible patch identities on one level while all four levels retain the
x-upper boundary branch.

Reuse the existing full-physics fresh and restart processes. Qualify serial and
sparse 1/2/4/8-rank field parity, independent serial continuation, and a
two-rank checkpoint restarted at four and eight ranks under changed ownership.
Keep schema 8 because the existing topology and fingerprint already represent
ordered parent/child branching without new compatibility state.

## 0.201.0 reproducible completion boundary

Freeze one recursive PeleC reference snapshot and define completion through
ordered capability workstreams with evidence-based exit gates. Keep project,
runtime, README, preset, mapping, and validation records under one automated
contract.

Reject MPI configurations that do not provide a usable Fortran 2008
`mpi_f08` module for the selected Fortran compiler. Replace host-associated
geometry callbacks with module procedures receiving explicit context, and
reject any produced ELF application or callback regression that requests an
executable stack.

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

## 0.202.0 three-dimensional coordinate core

Add the uniform 3D mesh contract and z-normal rotations around the existing
qualified x-normal Riemann implementations. Cover ideal-gas conserved and
primitive states, passive-multispecies closure, and general-EOS reacting
species flux placement with focused unit tests.

Do not expose a `pelef3d` application from this increment. The next increment
must add a conservative uniform-grid field update, CFL selection, periodic
boundaries, x/y/z dimensional-reduction tests, and a smooth convergence case
before moving to reactive transport or 3D AMR.

## 0.203.0 conservative periodic 3D Euler

Add a direct periodic x/y/z finite-volume divergence over the qualified
directional fluxes, a summed multidimensional CFL selector, and a
rollback-safe SSPRK2 update. Keep the first spatial boundary deliberately at
PCM so coordinate assembly, conservation, and time integration are isolated
from a future multidimensional high-order predictor.

Expose `pelef3d` with a strict constant-pressure diagonal entropy-wave input.
Require exact x/y/z whole-step reduction against the 1D PCM solver for both
Rusanov and PeleC-style Riemann paths, first-order smooth-wave convergence,
roundoff periodic conservation, and an independent deterministic CSV check.

Do not infer reacting-flow or production-3D completion from this gate. The
next regular-grid increment must extend the general-EOS reacting state,
chemistry splitting, and molecular transport before 3D AMR is introduced.

## 0.204.0 general-EOS multispecies 3D Euler

Extend the periodic 3D divergence and SSPRK2 transaction to the complete
runtime reactive state. Reuse the established general-EOS directional fluxes,
recover NASA7 temperature at each candidate stage, and integrate all Euler and
species components. Expose a separate `pelef_reactive_3d` application so the
constant-`gamma` baseline remains stable.

Qualify Rusanov, HLLC, and PeleC-style x/y/z whole-step reduction against an
independent 1D SSPRK2 reference, a seven-species entropy-wave refinement
study, roundoff all-component conservation, invalid-solver rollback, and a
full-H2O2-thermo public CSV contract.

Do not enable chemistry or molecular transport in this increment. Add the
cell-local reaction split next, then directional diffusive fluxes and their
timestep limit, before beginning 3D AMR.

## 0.205.0 cell-local chemistry splitting in 3D

Apply the established elementary explicit or full-H2O2 implicit
constant-volume reactor independently to every cell of one private 3D
candidate. Compose it symmetrically around the qualified 3D hydro step and
publish state and temperature only after both chemistry half-steps and hydro
succeed.

Qualify cell results against the existing one-dimensional reactor, uniform
whole-step reduction, H/O/N elemental conservation, nonzero reaction progress,
and exact outer rollback when hydro fails after the first chemistry half-step.
Expose a strict periodic Gaussian hotspot and independently check its dynamic
CSV schema, EOS identities, physicality, Euler/element conservation, and
coupled response.

Do not infer detailed-fuel, ignition-validation, or whole-PeleC parity from the
built-in H2/O2 gates. Add three-dimensional molecular transport and its stable
timestep next, before beginning 3D AMR.

## 0.206.0 regular-grid molecular transport in 3D

Extend the periodic general-EOS field with direct x/y/z diffusive faces. Reuse
the established mixture coefficient and species-flux kernels, but assemble the
complete three-dimensional Newtonian stress from all nine velocity-gradient
entries. Include Fourier conduction, mixture-averaged diffusion,
barodiffusion, correction velocity, species enthalpy transport, and one
six-face species-positivity bound per cell.

Advance transport with rollback-safe SSPRK2 and limit the public timestep by
the sum of the three inverse squared spacings. Compose optional chemistry,
transport, and hydro in one private `R-T-H-T-R` candidate. Require elementary
and full-H2O2 x/y/z reductions against the established 2D transport path,
second-order convergence of a genuinely three-dimensional viscous mode,
thermal/species smoothing and conservation, late-stage rollback, and a public
full-H2O2 transport/control comparison.

Keep this increment serial, periodic, uniform, and single level. Introduce 3D
AMR conservation and reflux next, followed by restart and distributed rank
parity before any 3D EB claim.

## 0.207.0 static two-level reactive 3D AMR hydro

Introduce one rectangular fine patch strictly inside a periodic coarse grid.
Use piecewise-constant prolongation for initialization, advance the provisional
coarse level with SSPRK2 while retaining its time-averaged x/y/z face fluxes,
and take `r` fine SSPRK2 substeps with coarse boundary states interpolated in
time. Average the `r^2` fine faces and `r` substeps at each interface, reflux
all six surrounding coarse faces, and volume-average the fine solution onto
covered parents.

Require roundoff composite conservation for every Euler and species component,
exact covered-parent synchronization, physical general-EOS recovery, uniform
state invariance, a measurably active reflux correction, and exact rollback
when a late face solve fails. Expose a strict full-H2O2 public entropy-wave
case with separate coarse/fine CSV output and independent contract checking.

Keep chemistry, molecular transport, high-order hydro/prolongation, dynamic
regridding, restart, MPI, boundary-touching refinement, and EB outside this
increment. Add restart identity next, then distributed ownership and
1/2/4/8-rank parity before composing AMR sources or diffusion.

## 0.208.0 static 3D AMR hierarchy checkpoint/restart

Add an explicit schema for the qualified two-level hierarchy. Persist both
conserved levels and temperatures only at synchronized root-step boundaries,
together with time, cumulative coarse steps, initial composite integrals, and
maximum reflux history. Fingerprint the full NASA7 species database and all
solver, mesh, patch, refinement, CFL, boundary, and enabled-physics settings
that can alter continuation.

Read into private candidates and publish nothing unless schema, extents,
metadata, physical temperature recovery, terminal record, and coarse/fine
average-down synchronization all pass. Reject checkpoint/output path
collisions before execution. Require missing, truncated, trailing, and
incompatible inputs to leave caller state and history exactly unchanged.

Run the public full-H2O2 entropy wave continuously and through a step-four
stop/restart boundary. Require byte-for-byte equality of both final level CSV
files. Keep this format serial, static, single-patch, hydro-only, and
schema-version-specific. Distributed ownership and 1/2/4/8-rank parity are
next; dynamic topology, AMR sources/diffusion, and 3D EB remain later work.

## 0.209.0 distributed-slab static 3D AMR hydro

Introduce a deterministic x-slab descriptor for both levels of the qualified
static hierarchy. Compute CFL rates, directional faces, and SSPRK2 cell
candidates only for owned planes; reconstruct complete candidates with fixed-
order collectives and reuse the exact serial reflux and average-down sequence.
Keep hierarchy memory replicated so distribution semantics can be qualified
before a sparse-field redesign.

Require fixed-size collective consensus over patch/layout metadata, solver and
timestep controls, and the complete NASA7 table. A root hierarchy broadcast
must reject inconsistent roots or extents before any variable-count payload
operation and must publish only a fully validated candidate.

Expose a namelist-driven MPI application using selected-root initialization,
checkpoint/restart, and CSV output. Require exact serial and 1/2/4/8-rank
coarse/fine outputs, plus a two-rank checkpoint resumed at four and eight
ranks. Defer sparse hierarchy storage, scalable I/O, AMR chemistry/transport,
high-order reconstruction, dynamic regridding, boundary-touching patches, and
3D EB.

## 0.210.0 regular and static-AMR characteristic PLM

Add a selectable frozen-composition characteristic-PLM spatial operator to
the periodic general-EOS 3D solver. Rotate y/z differences into the established
x-normal characteristic basis, apply MC or minmod limiting plus a common
physical-state scale, convert reconstructed faces through the NASA7 EOS, and
retain transactional SSPRK2 publication with time-averaged face fluxes.

Apply the same operator to both levels of the serial static hierarchy. Fill
two fine ghost layers from coarse endpoint states interpolated in time and
limited linearly in space at each fine-center offset. Feed the resulting
interface fluxes through the existing area/time average, reflux, and
average-down sequence. Extend the checkpoint fingerprint to reconstruction
and limiter controls and require exact continuous/restarted PLM outputs.

Require x/y/z symmetry, second-order smooth-wave convergence, roundoff
conservation, physical/transactional failures, and a measurable PLM accuracy
improvement over PCM on the public AMR entropy wave. Keep the distributed AMR
kernel PCM-only until its reconstruction-width halo contract is implemented;
defer CTU/PPM, sources/diffusion across levels, dynamic topology, physical
boundaries, and 3D EB.

## 0.211.0 static axis-plane 3D EB foundation

Introduce a standalone 3D EB metric type before porting redistribution or AMR
ownership. Build one analytic x-, y-, or z-normal plane with fluid on the
positive side. Retain cell volume fraction/centroid, all Cartesian face
apertures and tangential centroids, and EB area/centroid/normal data. Make the
cellwise discrete aperture-divergence/EB-normal identity a validator invariant
so orientation and units are fixed before general geometry is attempted.

Add a NASA7 slip-wall pressure flux and the corresponding aperture-weighted
finite-volume divergence. Build PCM directional Riemann fluxes on open faces,
use a zero-gradient exterior state at the physical domain, and recover active
cell temperature before publishing a forward-Euler candidate. Keep covered
cells untouched and make every error path transactional.

Require exact x/y/z uniform tangent-flow invariance, integrated pressure
force, no wall mass/energy/species leakage, and roundoff all-component
fluid-volume conservation for a nonuniform field. Explicitly reject an
interior face-aligned plane until boundary-metric ownership exists. Add
small-cell redistribution and its stable timestep next; only then broaden the
geometry and compose 3D EB transport, AMR/reflux, restart, and MPI ownership.

## 0.212.0 conservative planar 3D EB FluxRedist

Port the established first-order FluxRedist conservation rule to a
six-connected 3D neighborhood. Use stored face apertures to admit only active
neighbors, volume fractions for the neighborhood average and excess return,
and private candidates for transactional publication. Keep the operation
component-agnostic so Euler and every species density follow identical
arithmetic.

Compose the operator first with a caller-supplied residual and NASA7 recovery,
then with the qualified PCM open-face and slip-wall divergence. Preserve the
raw Euler routine as a control. Require analytical x/y/z one-neighbor values,
overlapping transverse-neighbor conservation, uniform residual identity, and
a full-grid CFL case where raw density is negative but redistribution is
physical.

Do not call this weighted StateRedist or broad PeleC EB parity. Implement the
weighted/higher-order state algorithm and a public EB timestep/application
contract next, before generalizing geometry or adding 3D EB transport,
AMR/reflux/regridding, restart, and MPI ownership.

## 0.213.0 zeroth-order weighted StateRedist for planar 3D EB

Keep FluxRedist intact and add a separately named provisional-state path.
For the qualified exact axis plane, select the one regular cell along the
aperture-difference normal. Raise every small neighborhood to a configurable
target volume, partition shared recipients by their neighborhood count, and
use the same weights for volume gather and result scatter so every conserved
component remains invariant in `sum(kappa U)`.

Apply the order-zero algorithm first to arbitrary provisional fields, then to
a caller-supplied reactive residual and NASA7 recovery. Refactor the integrated
hydro routes to consume one common face-flux and EB-divergence builder. Require
x/y/z closed-form weights, uniform identity, covered-state behavior,
conservation, positivity restoration, custom-target control, and complete
transaction rollback. Repeat the full-grid raw-negative stability control for
both FluxRedist and StateRedist.

Reject non-axis normals or missing regular receivers instead of implying
general geometry. Keep higher-order StateRedist, a public timestep/application,
oblique/curved geometry, characteristic EB reconstruction, chemistry,
transport, AMR/reflux/regridding, restart, and MPI ownership as subsequent
increments.

## 0.214.0 public stabilized planar 3D EB hydro

Add a full-grid CFL selector over active general-EOS cells and pair it only
with FluxRedist or StateRedist. Introduce a dedicated validated public input
record for one exact axis plane with a regular receiver; exclude the raw route
from public selection.

Install a serial application that initializes the frozen elementary-mixture
density sheet, advances repeated transactional forward-Euler steps, clips the
last timestep, computes fluid-volume diagnostics, and emits all geometry and
reactive fields in deterministic x-fastest CSV. Separate wall-normal pressure
impulse from the genuinely conserved fluid invariants.

Run both stabilization methods through the public executable and independently
reconstruct geometry, EOS/composition closure, integrals, response, and
method-specific signatures from their outputs. Defer characteristic/higher-
order EB evolution, chemistry/transport, general geometry, AMR/reflux,
restart, MPI ownership, and scalable output.

## 0.215.0 optional elementary chemistry in planar 3D EB

Extend the transactional 3D cell-local chemistry wrapper with an optional
active mask and route every z plane through the already qualified 2D reactor.
Build the EB mask strictly from noncovered cells and require exact covered
storage after each candidate stage.

Compose the reaction operator around the existing stabilized hydro step as
`R(dt/2)-H(dt)-R(dt/2)`. Keep the full split private until both chemistry
halves and FluxRedist or StateRedist hydro succeed. Preserve the disabled path
bitwise and reject empty chemistry, invalid tolerances, or an invalid hydro
selection transactionally.

Expose only the frozen elementary H2/O2/N2 mechanism in the public application
and add H/O/N fluid-volume diagnostics. Qualify one reacting StateRedist case
against a byte-identical inert `0.214.0` baseline, frozen response signatures,
and independent geometry/EOS/conservation checks. Keep the pre-reaction CFL
and first-order PCM restriction explicit; defer reaction-aware timestep
control, configurable/stiff chemistry, public FluxRedist chemistry,
transport, higher-order/general EB, AMR/reflux/regridding, restart, MPI, and
scalable output.

## 0.216.0 planar 3D EB molecular transport

Expose the existing regular-grid three-dimensional constitutive face law to a
dedicated EB transport module. Recover only active primitives, evaluate open
interior Cartesian faces, impose zero transport flux on physical-domain and
adiabatic slip/impermeable EB faces, and form the aperture/physical-volume
divergence. Keep covered storage outside EOS and timestep scans.

Apply a physical species-inventory outflow limiter before divergence, including
the matching species-enthalpy energy correction. Advance each half transport
interval with SSPRK2 and order-zero StateRedist at both Euler stages. Compose
reaction, transport, and PCM hydro as a private `R-T-H-T-R` transaction;
preserve the existing `R-H-R` result bitwise when transport is disabled.

Add default-off validated switches for viscosity, conduction, species
diffusion, and barodiffusion plus a transport CFL. Select the minimum active-
state hydro and transport bounds. Qualify x/y/z conservation, zero boundary
and wall transport, limiter activation, covered identity, rollback, and
disabled parity. Freeze one-step control/transport/coupled public outputs with
independent geometry, EOS, relative-invariant, element, response, schedule,
diagnostic, and hash checks.

Keep the boundary narrow: serial, single level, exact axis plane, PCM hydro,
StateRedist transport, elementary seven-species mixture, zero outer transport,
and adiabatic slip/impermeable EB wall. Defer general wall transport,
FluxRedist transport, reaction-aware stage timestep recomputation, general
geometry/centroid interpolation, higher-order EB hydro, AMR/reflux/regridding,
restart, MPI ownership, and external validation.

## 0.217.0 transactional planar 3D EB checkpoint/restart

Add a dedicated versioned checkpoint for the serial single-level planar 3D EB
application. Store the complete thermodynamic, elementary-kinetics, and
transport records; immutable geometry and numerical controls; the full state
and temperature; initial conserved/L1/element inventories; and cumulative
time-step and transport diagnostics. Treat output paths, checkpoint cadence,
and continuation limits as runtime policy while requiring the stored clock to
fit the requested continuation.

Make restart publication transactional: parse and validate private candidates,
including the terminal marker and active-cell EOS temperature, before changing
any caller-owned state or metadata. Reject truncated, corrupted, or
physics-incompatible files with bitwise rollback.

Qualify a coupled chemistry-plus-transport run split after its first committed
step. Require a real continuation of at least one further step, exact final CSV
parity with the uninterrupted run, identical cumulative diagnostics, and
explicit checkpoint magic/schema/end-marker checks. Keep formatted serial I/O,
one exact axis plane, and the existing frozen elementary mixture; defer MPI or
AMR restart, scalable I/O, general geometry, and schema conversion.

## 0.218.0 pinned Cantera YAML ingestion

Pin the exact Cantera H2/O2 YAML source in the repository and require the
ideal-gas `ohmech` phase explicitly. Use Cantera's parsed API values while
retaining source hash, metadata/API versions, source and target units, element
composition, reaction indices, duplicate flags, and named-collider semantics.
Reject any phase or record family that cannot be represented faithfully by the
current Fortran kernels.

Extend the normalized JSON generator to emit NASA7 and transport loaders in
the same module as the existing reaction and Jacobian kernels. Replace only
the full-H2/O2 hand-maintained thermo and transport tables; keep the elementary
subset and its byte-clean generated source unchanged.

Require helper/unit tests with and without Cantera, byte-exact pinned
YAML-to-JSON and JSON-to-Fortran regeneration, exact reaction-family and
provenance checks, the existing full-H2/O2 CFD path, and live Cantera 0D and
1D/2D parity. Keep this a bounded build-time ingestion route; defer runtime
loading, CHEMKIN, unsupported thermo/rate families, arbitrary detailed fuels,
and production stiff integration.

## 0.219.0 selected normalized mechanism build

Add one optional CMake entry point that reads the names declared by a
normalized JSON bundle, regenerates its Fortran source in the build tree, and
compiles it separately from the fixed application library. Configure a probe
against those names so a selected bundle must execute thermo, kinetics,
transport, rate, and Jacobian interfaces rather than merely produce text.

Accept an optional source YAML path and verify its recorded basename and
SHA-256 both during configuration and immediately before regeneration. Keep
the default production target set unchanged and install the probe only in a
selected-bundle configuration.

Qualify the reusable CMake path with a small compiled fixture and the pinned
full-H2/O2 bundle in an independent tests-disabled build. Strengthen the
existing reduced reactor Jacobian with a directional finite-difference test as
preparation for the next stiff-integrator increment. Defer application-level
dispatch, runtime parsing, expanded mechanism forms/sizes, and CVODE itself.

## 0.220.0 selected normalized mechanism 0D reactor

Configure and install a second selected-bundle executable that loads the
generated NASA7 and reaction tables and advances the existing adaptive
implicit constant-volume reactor. Isolate the generic runtime modules from
the fixed generated mechanisms so a selected bundle may safely reuse their
module names.

Accept initial mole fractions by exact species name, map them into bundle
order, and reject nonfinite controls, malformed or duplicate compositions,
unknown names, invalid reactions, molecular-mass imbalance, and temperatures
outside the common NASA7 interval before creating output. Reject portable
lexical input/output aliases. Bound the initial step and output cadence against
the adaptive minimum while permitting an unreduced scheduler fragment to land
exactly on the final time. Extend generator validation to reaction element and
molecular-mass balance. Require the CMake helper to reject generator-only
legacy kinetics bundles unless schema-1 thermo/transport records and explicit
interfaces are complete.

Qualify the public path with a reactive two-species fixture, independent
conservation and activity checks, repeated-run byte identity, explicit startup
failure diagnostics with no CSV artifact, and byte-exact selected/fixed
full-H2/O2 trajectory parity. Add a below-minimum final-fragment regression.
Install and stack-audit the selected executable.
Keep runtime parsing and dispatch, fixed CFD application integration,
CVODE/SUNDIALS parity, mechanisms above 32 species, new rate families, and
detailed-fuel validation for later increments.

## 0.221.0 optional official-Fortran CVODE selected reactor

Add a default-off SUNDIALS boundary isolated from the native selected runtime.
Require exactly the pinned 7.2.0 official static Fortran targets for CVODE,
serial N_Vector, dense SUNMatrix, and the dense linear solver; reject C-only,
shared-only, different-version, and non-double-precision prefixes. Keep one
selected application and dispatch `integrator = "implicit"` or `"cvode"`
without silent fallback.

Advance `N-1` independent mass fractions with BDF and reconstruct the largest
initial species from closure. Reuse the generic fixed-density,
fixed-internal-energy RHS and generalize the energy-constrained reduced
Jacobian to an arbitrary dependent species. Keep the generated mechanism
kernel limited to reported production rates. Make caller publication
transactional, quarantine a failed solver, and enforce a distinct cumulative
CVODE internal-step budget across output calls.

Qualify the two-species closure map, forced and cumulative-limit rollback,
deterministic output, below-minimum final scheduling, pinned full-H2/O2
conservation, Cantera trajectory/rate comparison, and a tests-disabled static
install. Keep the native backend the default and defer multiple live contexts,
thread safety, sparse solvers, CFD dispatch, detailed fuels, and performance.

## 0.222.0 serial CVODE multi-context ownership

Replace the singleton API with a context-first opaque handle. Store only a
private slot and generation in that handle. Keep every noninteroperable
mechanism object and SUNDIALS resource in a bounded private registry, and pass
only a stable C-interoperable slot/generation token through
`FCVodeSetUserData`. Reject null, inactive, or generation-stale tokens before
callbacks touch Fortran state.

Permit up to 64 live contexts and arbitrary sequential interleaving. Keep
clocks, solver statistics, cumulative step budgets, temperature guesses, and
failure quarantine strictly local. Release resources in dependency order and
make same-handle finalization idempotent. Treat intrinsic handle copies as
aliases: the first valid finalization owns cleanup, and any surviving copy must
fail stale-generation validation rather than release a reused slot. Reject
capacity exhaustion without assigning the requested handle. If a vendor
destroy routine returns an error, continue dependency-ordered cleanup,
invalidate the handle, and report the first failing routine without promising
an unsafe retry.

Require exact trajectory and solver-statistics identity between standalone
and interleaved states with different density, temperature, and closure
species. Force one context to fail while its peer continues, finalize peers in
arbitrary order, exercise stale-slot reuse, reverse cleanup of all 64 slots,
and retain every `0.221.0` fixture, full-H2/O2, Cantera, native, and install
gate. Do not claim threaded/reentrant use, sparse integration, performance,
runtime mechanism loading, CFD dispatch, or broader chemistry validation.

## 0.223.0 configure-time selected regular 1D reacting flow

Isolate the regular 1D reacting-flow dependencies from `pelef_core` so a
selected generated mechanism can legally reuse the committed full-H2/O2
module and symbol names. Introduce one mechanism-independent primitive gas
transport type/builder and make the fixed database layer reexport it. Keep the
ordinary fixed application source-compatible and require its parser to reject
`chemistry_model = "selected"` unless a configured selected application opts
in explicitly.

Configure and install `pelef_reactive_1d_selected` alongside the selected
probe and 0D reactor. Load the selected NASA7, reaction, and transport arrays;
validate their structure, molecular-mass balance, finite common temperature
range, and requested hotspot range before output creation. Resolve initial
composition by exact nonduplicate species name and pass the normalized
bundle-order composition into the existing regular 1D initializer and native
chemistry/transport evolution.

Qualify an independently ordered reacting two-species fixture for activity,
closure, conservation, uniformity, deterministic output, and negative startup
behavior. Select the pinned ten-species/29-reaction H2/O2 bundle and require
byte-identical CSV output against the fixed regular-1D executable. Retain the
full eight-configuration regression matrix and perform a clean tests-disabled
build/install/ELF/dependency audit.

Reject AMR and checkpoint/restart controls in this application. Keep runtime
mechanism loading, CVODE-backed CFD chemistry, selected 2D/3D/AMR/EB/MPI
dispatch, mechanisms above 32 species, thread safety and performance,
detailed-fuel validation, and external PeleC/experiment parity for later
increments.

## 0.224.0 configure-time selected regular 3D reacting flow

Move exact-name composition resolution and selected-mechanism structural,
balance, transport-order, and common-temperature-range checks into shared
startup modules used by every selected CFD front end. Keep the ordinary fixed
3D parser closed to `thermo_model = "selected"` unless a configured selected
application opts in explicitly.

Extract the regular 3D mesh, initialization, CFL/transport timestep,
transactional `R-T-H-T-R` loop, invariant checks, CSV publication, and
diagnostics into one mechanism-independent application driver. Configure and
install `pelef_reactive_3d_selected` with a private link graph containing the
selected generated bundle and generic thermo, kinetics, transport, and regular
3D modules, but neither committed generated mechanism.

Qualify an independently ordered two-species 4-by-4-by-4 periodic fixture for
active chemistry, positivity, closure, conservation, spatial uniformity, and
repeat-run byte identity. Select the pinned ten-species/29-reaction full-H2/O2
bundle and require byte-identical transport-hotspot CSV against the fixed 3D
front end, while retaining the independent nonzero transport-response check.
Reject unknown composition names and selected input sent to the fixed front
end before output creation. Add a short nonuniform full-H2/O2 chemistry case
using characteristic PLM and require byte-identical selected/fixed output so
the selected full reaction family and higher-order dispatch are both live.

Retain the complete eight-configuration regression matrix and perform a fresh
tests-disabled selected Release build, install, ELF-stack, dependency, and
build/install identity audit. Keep physical boundaries, selected 2D/AMR/EB/MPI,
runtime mechanism loading, CVODE-backed CFD chemistry, mechanisms above 32
species, thread safety and performance, detailed-fuel validation, and external
PeleC/experiment parity for later increments.

## 0.225.0 configure-time selected regular 2D reacting flow

Extract the fixed regular 2D simulation, deterministic CSV output, invariant
calculation, extrema, and diagnostics into one mechanism-independent
application driver. Keep the fixed executable responsible only for fixed
model loading and fixed composition resolution. Configure and install
`pelef_reactive_2d_selected` from the normalized bundle's declared interfaces
and SHA-256.

Build a private selected 2D runtime on the selected 1D core. Compile only the
regular 2D mesh, configuration, physical-boundary, CTU, native chemistry,
molecular-transport, and application modules, and keep `pelef_core`, fixed
database loaders, and committed generated mechanisms out of the link graph.
Make selected 3D reuse this target instead of recompiling the same 2D modules.

Extend the regular 2D configuration with bounded exact-name composition
fields and 32-species wall-flux storage. Preserve fixed parser behavior unless
the configured selected application explicitly enables `selected`. Validate
the actual selected species count after loading so every unused wall-flux slot
must be zero. Thread the resolved bundle-order composition through
initialization, composition-wave, boundary, and simulation paths without
changing existing fixed callers.

Qualify an independently ordered two-species fixture for active chemistry,
positive state, closure, uniformity, and repeat-byte identity. Require
byte-identical fixed/selected pinned full-H2/O2 output for uniform implicit
chemistry, a nonuniform characteristic-PLM/CTU active-chemistry step, and
prescribed-species-wall molecular transport. Reject unknown/duplicate names,
nonfinite or nonzero tails, unsafe composition-wave endpoints, selected input
sent to the fixed executable, overlapping hotspots outside the common NASA7
range, and isothermal ghost temperatures outside that range before output.

Retain the complete eight-configuration matrix and perform a clean
tests-disabled selected Release build/install/ELF/dependency audit. Keep
selected AMR, EB, MPI, runtime mechanism loading, CVODE inside CFD, mechanisms
above 32 species, thread safety and performance, detailed-fuel validation, and
external PeleC/experiment parity for later increments.

## 0.230.0 configure-time selected serial reactive EB 3D flow

Extract the fixed serial reactive EB 3D application's exact planar geometry,
adaptive CFL loop, optional `R-T-H-T-R` advancement, redistribution,
deterministic CSV publication, invariants, and diagnostics into one
mechanism-independent driver. Preserve the fixed namelist and freeze the
coupled seven-species output before extraction.

Build a private selected EB 3D runtime on the selected regular-3D core.
Compile only the generic geometry, CFL, wall-flux, redistribution, hydro,
transport, driver, checkpoint-support dependency, and shared application
modules; keep `pelef_core`, fixed loaders, and committed generated mechanisms
out of its link graph. Thread optional bundle-order composition through
regular/cut initialization and the generated chemistry policy through both
masked chemistry half-steps.

Qualify an independently ordered H2/H exact-plane fixture for active
chemistry, geometry, physicality, closure, exact covered storage, and repeated
byte identity. Require a complete test-only seven-species selected bundle to
match fixed output byte-for-byte, and exercise the pinned full-H2/O2 bundle
with implicit chemistry and all qualified molecular transport. Reject invalid
composition, temperature, path, persistence, element-family, and fixed-parser
inputs before output.

Retain the complete eight-configuration matrix and perform a fresh
tests-disabled Release build/install/ELF/dependency audit with all 38
executables. Keep selected checkpoint/restart, general element diagnostics
and EB geometry, EB AMR/MPI 3D, runtime mechanism loading, CVODE inside CFD,
mechanisms above 32 species, thread safety and scaling, detailed-fuel
validation, and external PeleC/experiment parity for later increments.

## 0.231.0 selected serial reactive EB 3D checkpoint/restart

Preserve fixed planar EB checkpoint schema 1 byte-for-byte. Add an exclusive
schema 2 for selected calls and place a selected-context record before the
existing complete species/reaction/transport/configuration/geometry/state
body. Store the configure-time bundle SHA-256, generated chemistry-integrator
policy, species count, and normalized bundle-order mole fractions. Require
all selected context fields together; do not permit selected-to-fixed schema
downgrade or cross-runtime schema reads.

Connect the already supplied selected SHA, composition, and policy through the
shared application checkpoint calls. Enable checkpoint cadence,
stop-after-write, and restart in the selected front end only after that context
is present. Reject input/checkpoint and input/restart lexical aliases before a
write or read can touch the input.

Extend the transactional checkpoint unit gate for selected round trip,
cross-schema rejection, changed SHA/policy/composition, malformed or truncated
context, nonfinite/negative composition, and complete rollback. Add a
separate-process selected elementary uninterrupted/checkpoint-stop/restart
chain with byte-exact final CSV and diagnostic parity plus an application-level
composition-mismatch rejection. Recheck the frozen fixed schema-1 hashes.

Retain the complete regression/build/install matrix. Keep checkpoint payload
digests, crash-atomic replacement, schema conversion, general geometry,
selected AMR/MPI EB 3D restart, runtime mechanism loading, CFD CVODE, detailed
fuels, thread safety, performance, and external validation for later
increments.

## 0.229.0 configure-time selected serial EB AMR 2D flow

Extract the fixed serial EB AMR application's hierarchy setup, subcycled
`R-T-H-T-R` loop, reflux, average-down, deterministic level output,
invariants, and diagnostics into one mechanism-independent driver. Preserve
the fixed namelist and numerical path. Let fixed and selected front ends load
their own models and compositions.

Build a private selected EB AMR runtime on the selected single-level EB core.
Compile only the generic hierarchy, reflux, regrid, transport, driver, and
shared application modules; keep `pelef_core`, fixed loaders, and committed
generated mechanisms out of its link graph. Thread optional bundle-order
composition through coarse/fine and boundary initialization and the generated
chemistry policy through both half-steps on both levels.

Qualify a static two-level independently ordered fixture for geometry,
physicality, closure, and exact repeat output. Require selected implicit
full-H2/O2 coarse and fine results to be byte-identical to the fixed
executable while chemistry, molecular transport, characteristic PLM, and an
isothermal moving no-slip embedded wall are active. Reject dynamic,
three-level, multipatch, checkpoint/restart, path alias, unknown-species, and
out-of-range-temperature inputs before output creation.

Retain the complete eight-configuration matrix and perform a sanitized,
tests-disabled selected MPI Release build/install audit with all 37
executables. Keep selected EB 3D, dynamic/multipatch EB AMR, MPI AMR/EB,
checkpoint/restart, runtime mechanism loading, CVODE inside CFD, mechanisms
above 32 species, thread safety and scaling, detailed-fuel validation, and
external PeleC/experiment parity for later increments.

## 0.226.0 configure-time selected serial AMR 1D reacting flow

Extract the fixed serial reactive 1D AMR mode dispatch, simulation, composite
CSV output, conservation calculation, topology reporting, and diagnostics into
one mechanism-independent application driver. Keep the fixed executable
responsible for committed-model loading and fixed composition resolution.
Configure and install `pelef_amr_reactive_1d_selected` from the normalized
bundle's declared interfaces and SHA-256.

Build a private selected AMR runtime on the selected regular-1D core. Compile
only the hierarchy, multipatch, regrid, reacting AMR, and application modules;
keep `pelef_core`, fixed database loaders, and committed generated mechanisms
out of the link graph. Thread an optional bundle-order root composition through
two-level, arbitrary-depth, and multipatch initialization without changing
fixed call sites.

Qualify an independently ordered two-species fixture for active chemistry,
transport, refinement, positive state, closure, composite coverage, and
repeat-byte identity. Require byte-identical fixed/selected full-H2/O2 output
for two-level active chemistry/transport, three-level characteristic PPM, and
three-patch molecular transport. Reject unknown species, non-AMR input,
out-of-range entropy-wave temperature extrema, and selected input sent to the
fixed executable before output.

Retain the complete eight-configuration matrix and perform a clean
tests-disabled selected Release build/install/ELF/dependency audit. Keep
selected AMR checkpoint/restart, EB, MPI, runtime mechanism loading, CVODE
inside CFD, mechanisms above 32 species, thread safety and performance,
detailed-fuel validation, and external PeleC/experiment parity for later
increments.

## 0.227.0 configure-time selected serial reactive EB 2D flow

Extract the fixed single-level reactive 2D EB geometry/boundary construction,
simulation, deterministic CSV publication, volume-weighted invariant checks,
and diagnostics into one mechanism-independent application driver. Keep fixed
and selected front ends responsible for their own model loading, exact-name
composition resolution, and startup temperature validation.

Build a private selected EB runtime on the selected regular-2D core. Compile
only the generic EB geometry, reconstruction, hydro, transport, driver, and
application modules; keep `pelef_core`, fixed loaders, and committed generated
mechanisms out of the link graph. Thread optional bundle-order composition
through both regular initialization and configured physical-boundary states,
and thread the generated chemistry-integrator policy through both masked
Strang chemistry stages without changing fixed call sites.

Qualify an independently ordered two-species plane fixture for active
chemistry and transport, positive active state, mass-fraction closure,
regular/cut/covered geometry counts, and repeat-byte identity. Require
byte-identical fixed/selected full-H2/O2 output for active plane chemistry and
for a characteristic-PLM circle with molecular transport and an isothermal
moving no-slip embedded wall. Reject unknown species, an out-of-range
embedded-wall temperature, a lexical input/output alias, and selected input
sent to the fixed executable before output creation.

Retain the complete eight-configuration matrix and perform a clean
tests-disabled selected Release build/install/ELF/dependency audit with all
seven selected executables. Keep selected EB AMR/3D, checkpoint/restart, MPI,
runtime mechanism loading, CVODE inside CFD, mechanisms above 32 species,
thread safety and performance, detailed-fuel validation, and external
PeleC/experiment parity for later increments.

## 0.228.0 configure-time selected regular MPI 1D flow

Extract the fixed coupled MPI verification application's decomposition,
initialization, adaptive `R-T-H-T-R` loop, collective diagnostics, gather, and
CSV publication into one mechanism-independent driver. Keep fixed and selected
front ends responsible for their own model loading. Preserve the fixed
ten-species initialization, command line, and CSV byte-for-byte.

Build a private selected MPI runtime on the selected regular-1D core and
`MPI::MPI_Fortran`. Compile only the generic MPI domain, reactive transport,
reactive advance, and shared application modules; keep `pelef_core`, fixed
loaders, and committed generated mechanisms out of its link graph. Thread the
optional generated chemistry-integrator policy through both distributed
chemistry half-steps and adaptive retries. Require communicator-wide policy
consensus and collective rollback on invalid or rank-disagreeing input.

Qualify an independently ordered explicit two-species fixture for activity,
physicality, closure, conservation, and exact 1/2/4-rank output. Require
selected implicit full-H2/O2 output to be byte-identical across 1/2/4 ranks,
to the fixed executable, and to the frozen pre-refactor one-rank output. Run
invalid/rank-disagreeing policy rollback at all three rank counts and
generalize the comparator to dynamic species columns. Keep the selected MPI
initializer bounded to the exact ordered H2/H fixture and pinned ten-species
full-H2/O2 profiles.

Retain the complete eight-configuration matrix and perform a sanitized,
tests-disabled selected MPI Release build/install/ELF/dependency audit with
all 36 executables. Keep general input-driven MPI flow, selected MPI AMR/EB,
checkpoint/restart, runtime mechanism loading, CVODE inside CFD, mechanisms
above 32 species, thread safety and scaling, detailed-fuel validation, and
external PeleC/experiment parity for later increments.

## 0.232.0 selected static two-level reactive EB AMR 2D checkpoint/restart

Preserve fixed two-level checkpoint schema 3 and its checkpoint-stop bytes.
Add exclusive selected schema 4 with an all-or-none context record containing
the configure-time bundle SHA-256, generated chemistry-integrator policy,
species count, and normalized bundle-order mole fractions. Reject partial
context, schema downgrade, and both cross-schema directions transactionally.

Thread the selected context through the shared static two-level application
without enabling dynamic, three-level, or multipatch persistence. Propagate
specific read/write failure context to the selected front end. Permit
checkpoint cadence, intentional stop, and restart, while rejecting lexical
aliases among the input, coarse output, fine output, checkpoint, and restart
paths before any file is touched.

Extend the checkpoint unit gate for selected round trip, incomplete write and
read context, fixed/selected cross-schema rejection, changed SHA/policy/
composition, malformed marker, truncated context, invalid composition,
rollback, and non-destructive failed writes. Add a separate-process pinned
full-H2/O2 fixed reference, selected reference, first-step checkpoint-stop,
restart, and composition-mismatch chain. Require exact fixed/selected and
uninterrupted/restarted coarse and fine CSV, schema-4 context validation, and
frozen reference, stopped, and checkpoint SHA-256 values. Also freeze the
unchanged fixed schema-3 checkpoint and stopped-CSV hashes.

Retain the complete regression/build/install matrix. Keep selected dynamic,
three-level, multipatch, and MPI EB AMR persistence, runtime mechanism
loading, crash-atomic and scalable I/O, payload authentication, CFD CVODE,
detailed fuels, thread safety, performance, and external validation for later
increments.

## 0.233.0 selected dynamic two-level reactive EB AMR 2D checkpoint/restart

Preserve fixed schema 3 and selected static schema 4 byte-for-byte. Add
exclusive selected dynamic schema 5 with the same all-or-none bundle SHA-256,
generated integrator, species-count, and normalized composition context plus
the initial composite-integral baseline.
Select the expected schema from the validated fixed/selected and
static/dynamic call contract, and reject both static/dynamic cross-schema
directions without fallback.

Permit `dynamic_regridding` in the selected front end only for the established
serial exactly-two-level, single-fine-patch application. Continue to reject
three-level, multipatch, and dynamic-parent configurations. Use the existing
body records to persist the actual post-regrid patch, regrid cadence and count,
clock, diagnostics, geometry, and coarse/fine fields transactionally. Restore
the schema-5 baseline before advancing so cumulative conservation diagnostics
remain identical to uninterrupted execution.

Extend the unit gate with schema-5 baseline round trip and corruption,
static/schema-5 and dynamic/schema-4 rejection, rollback, and unchanged
schema-4 coverage. Add
separate-process full-H2/O2 fixed and selected uninterrupted references, a
topology-changing first-step checkpoint-stop, independent restart, changed
composition rejection, exact fixed/selected and uninterrupted/restarted CSV
comparison, and frozen schema-5/checkpoint-stop hashes.

Retain the complete regression/build/install matrix and installed dynamic
restart smoke. Keep selected three-level, multipatch, dynamic-parent, and MPI
EB AMR persistence, runtime loading, crash-atomic/scalable I/O, payload
authentication, CFD CVODE, detailed fuels, thread safety, performance, and
external validation for later increments.

## 0.237.0 selected sparse MPI reactive AMR 1D checkpoint/restart

Preserve fixed patch-tree magic/schema 1 byte-for-byte. Add selected schema 2
under the same magic with the complete generated bundle SHA-256, integrator,
species count, normalized bundle-order composition, and original composite-
integral baseline before the established rank-neutral payload. Derive the
expected schema from the call contract and reject cross-schema input.

Build the selected target only with MPI and a complete selected bundle. Share
the application lifecycle with the fixed front end, require communicator-wide
agreement on selected context before state creation, initialize from the
selected composition, and forward the generated integrator through every
sparse chemistry half-step. Keep rank count and owner assignments outside the
checkpoint so restart rebuilds ownership for the current communicator.

Extend units with context/baseline round trip, fixed/selected cross-read
rejection, mismatch/truncation/trailing-content rollback, invalid-write non-
destruction, and one-/two-/four-rank integrator-policy propagation. Add a
separate-process pinned full-H2/O2 fixed reference, selected one-/two-/four-
rank references, one-rank checkpoint-stop, two-/four-rank restart, five
mutated-context failures, six physical payload/geometry failures including a
near-floor negative species density, rank-consensus failure, and startup path-
alias gate. Require byte-exact fixed/selected, rank-count, and restart outputs.

Retain all frozen restart gates, the complete configuration matrix, and a
clean installed Release smoke. Keep selected MPI EB AMR, runtime loading,
crash-atomic/scalable I/O, payload authentication, CFD CVODE, detailed fuels,
thread safety, performance, and external validation for later increments.

## 0.234.0 selected static three-level reactive EB AMR 2D checkpoint/restart

Preserve fixed static-three-level magic/schema 3 and fixed dynamic-three-level
magic/schema 4 byte-for-byte. Under the static magic, add selected schema 4
with the complete generated bundle SHA-256, integrator, species-count, and
normalized composition context plus the original composite-integral baseline.
Reject fixed/selected cross-schema reads without fallback.

Admit selected three-level execution only when regridding is static and the
hierarchy contains one root, one middle patch, and one finest patch. Forward
bundle-order composition through all three initial and boundary states and
forward the generated integrator through every chemistry half-step. Keep
selected dynamic-three-level, multipatch, and dynamic-parent startup
rejection. Include the finest output in every input/output/checkpoint/restart
lexical-alias check.

Extend the driver unit gate with selected three-level context/baseline round
trip, fixed/selected cross-schema rejection, bundle/integrator/composition
mismatch, corrupted baseline rejection, and empty/default failure outputs. Add a
separate-process pinned full-H2/O2 fixed reference, selected reference,
first-step checkpoint-stop, independent restart, and changed-composition
chain. Require exact fixed/selected and uninterrupted/restarted bytes on root,
middle, and finest levels; identical cumulative conservation diagnostics;
physicality and species closure; and frozen reference, stopped, and checkpoint
SHA-256 values.

Retain all frozen fixed static/dynamic and selected two-level restart gates,
the complete configuration matrix, and a clean installed Release smoke. Keep
selected dynamic-three-level, multipatch, dynamic-parent, and MPI EB AMR
persistence, runtime loading, crash-atomic/scalable I/O, payload
authentication, CFD CVODE, detailed fuels, thread safety, performance, and
external validation for later increments.

## 0.235.0 selected dynamic three-level reactive EB AMR 2D checkpoint/restart

Preserve fixed dynamic-three-level magic/schema 4 byte-for-byte. Add selected
schema 5 under the same dynamic magic with the complete generated bundle
SHA-256, integrator, species count, normalized composition, and original
composite-integral baseline before the unchanged dynamic three-level body.
Derive the required schema from the fixed/selected call contract and reject
schema or static/dynamic magic mismatch without fallback.

Admit the established serial exactly-three-level, one-patch-per-level selected
path with either finest-only or dynamic-parent regridding. Persist the actual
middle and finest patches and the dynamic policy, cadence, clock, regrid
history, diagnostics, and all three fields transactionally. Restore the
persisted patches and original baseline before continuation. Keep multipatch
selected persistence rejected.

Extend the driver unit gate with schema-5 round trip, baseline closure and
corruption, fixed/selected cross-schema and static/dynamic cross-magic
rejection, corrupted dynamic controls/patch/field/terminal records, invalid-
write non-destruction, and empty/default failed-read targets. Add a separate-
process pinned full-H2/O2 fixed reference, selected reference, dynamic-parent
topology-changing checkpoint-stop, independent restart, and changed-
composition chain. Require exact fixed/selected and uninterrupted/restarted
bytes at all three levels; identical cumulative conservation diagnostics;
physicality, chemistry activity, and species closure; and frozen reference,
stopped, and checkpoint SHA-256 values.

Retain all frozen fixed and earlier selected restart gates, the complete
configuration matrix, and a clean installed Release smoke. Keep selected
multipatch and MPI EB AMR persistence, runtime loading, crash-atomic/scalable
I/O, payload authentication, CFD CVODE, detailed fuels, thread safety,
performance, and external validation for later increments.

## 0.236.0 selected dynamic multipatch reactive EB AMR 2D checkpoint/restart

Preserve fixed patch-set magic/schema 3 byte-for-byte. Add selected schema 4
under the same magic with the complete generated bundle SHA-256, integrator,
species count, normalized composition, and original composite-integral
baseline before the unchanged patch-set body. Derive the required schema from
the fixed/selected call contract and reject cross-schema input without
fallback.

Admit the established serial dynamic exactly-two-level patch-set path. Persist
the complete committed child-patch set, regrid cadence and history, clock,
diagnostics, geometry, root field, and every child field transactionally.
Restore the persisted patch set and original baseline before continuation.
Forward bundle-order composition to root and child initialization and boundary
states, and forward the generated integrator through both chemistry half-steps
on every level. Keep selected MPI EB AMR persistence separate.

Extend the driver unit gate with schema-4 context/baseline round trip,
fixed/selected cross-schema rejection, composition mismatch, corrupted end
marker, trailing-content rejection, baseline corruption, invalid-write non-
destruction, empty/default failed-read targets, and direct integrator
propagation. Before mechanism loading, make the selected front end reject a
lexical alias between any input/output/checkpoint/restart path and every
possible derived `_patchNNNN` child output.

Add a separate-process pinned full-H2/O2 fixed reference, selected reference,
first-step checkpoint-stop, independent restart, and changed-composition chain
with two disjoint fine patches. Require exact fixed/selected and uninterrupted/
restarted bytes for the root and both children; identical cumulative
conservation diagnostics; physicality, chemistry activity, and species
closure; complete checkpoint parsing through end-of-stream; and frozen
reference, stopped, schema-3, and schema-4 checkpoint SHA-256 values.

Retain all frozen fixed and earlier selected restart gates, the complete
configuration matrix, and a clean installed Release smoke. Keep selected MPI
EB AMR persistence, runtime loading, crash-atomic/scalable I/O, payload
authentication, CFD CVODE, detailed fuels, thread safety, performance, and
external validation for later increments.

## 0.238.0 selected sparse MPI reactive EB patch-tree 2D execution

Extract the established fixed sparse MPI reactive EB patch-tree application
into a mechanism-independent driver. Keep fixed positional call compatibility
and schema-8 behavior. Add a configure-time selected front end that supplies
generated thermo, kinetics, transport, bundle SHA-256, bundle-order
composition, and integrator through the same lifecycle.

Extend sparse MPI EB chemistry, full-physics, and to-time APIs with one
trailing optional integrator. Resolve and consensus-check the policy before
candidate mutation and forward it through both chemistry half-steps. Validate
the complete selected context communicator-wide before owner-only root state
or boundary initialization.

Add a four-level dynamically regridded full-H2/O2 fixed/selected parity case,
one-/two-/four-rank exact-output gates, direct no-mutation policy tests, and a
rank-context mismatch process. Reject selected checkpoint cadence, stop,
checkpoint/restart paths, model mismatch, and path aliasing before output.
Freeze the fixed GNU Debug/Release schema-8 hashes and compare the Release
bytes directly with the previous milestone executable.

Compile and execute the same selected MPI EB front end from the two-species
fixture bundle as a separate explicit-integrator genericity gate. Validate its
four-level sparse CSV independently rather than relying only on the pinned
full-H2/O2 implicit path.

Retain the complete configuration matrix and clean installed Release/ELF/MPI
audit. Defer selected sparse MPI EB persistence to a separate context-bound
schema increment; also defer fixed-depth MPI modes, runtime loading, scalable
I/O, performance, detailed-fuel validation, and external PeleC parity.

## 0.239.0 selected sparse MPI reactive EB patch-tree 2D restart

Preserve fixed magic/schema 8 byte-for-byte. Add exclusive selected schema 9
to the shared serial patch-tree checkpoint body. Bind the complete generated
bundle SHA-256, integrator, normalized bundle-order composition, existing
schema-8 fingerprint, original composite-integral baseline, clock, operator
counters, regrid history, topology, and every patch field. Require complete
selected metadata on write and strict end-of-stream on read.

Validate selected context independently in the sparse MPI I/O boundary before
gather or read. Store no communicator size or owner assignment. Continue to
gather only to the I/O root, then on restart broadcast validated rank-neutral
topology, construct distribution for the active communicator, and scatter
fields directly to current owners. Explicitly branch fixed and selected I/O
in the shared application and admit selected persistence in the front end only
after all input/output/checkpoint/restart aliases have been rejected.

Add schema-9 unit round trip, cross-schema and context mismatch, invalid-
baseline write non-replacement, fixed schema-8 freeze, one-/two-/four-rank
uninterrupted references, one-rank checkpoint-stop, independent two-/four-
rank restart, exact-output checker, 16 payload/context corruption gates, and
six lexical-alias gates. Directly exercise rank-dependent selected metadata
presence on both MPI read and write boundaries and require collective
rejection before transfer or publication. Pin GNU Debug and Release checkpoint/final/stopped
SHA-256 values, rerun the full configuration matrix, and perform a clean
tests-disabled installed Release restart smoke.

Do not claim payload authentication, crash-atomic or scalable I/O, schema
migration, runtime mechanism loading, fixed-depth MPI modes, CFD CVODE,
performance/thread qualification, detailed-fuel physics, or external PeleC
field parity from these gates.

## 0.240.0 rank-local sparse static 3D AMR and distributed PLM

Replace the public replicated 3D AMR hierarchy with uneven rank-local coarse
x slabs and fine slabs aligned to parent ownership. Permit ranks with zero
fine cells. Propagate two periodic coarse ghost planes, exchange fine internal
x ghosts only between active owners, and reconstruct exterior fine ghosts
from coarse start/end halos at the fine-stage time.

Refactor the regular and fine-patch characteristic-PLM interfaces so they can
evaluate only an owned slab without allocating or reconstructing a global
field. Keep PCM and PLM SSPRK2 state, interface fluxes, reflux, average-down,
and temperature recovery private until communicator-wide acceptance. Compare
ownership arrays, patch layout, NASA7 records, floating controls, and solver/
reconstruction/limiter text before data-dependent communication.

Retain schema-2 bytes and rank-neutral persistence by gathering a complete
hierarchy only on the selected I/O root and scattering it after initialization
or restart. Require exact serial versus one-/two-/four-/eight-rank PCM and PLM
output, direct sparse halo and empty-fine-rank tests, collective metadata
rollback, and a two-rank checkpoint resumed independently on four and eight
ranks.

Do not describe root-formatted gather/scatter as scalable I/O. Keep dynamic
topology, AMR chemistry/transport, physical boundaries, CTU/PPM, 3D EB AMR,
performance scaling, and external PeleC field parity for later increments.

## 0.241.0 fixed chemistry in static 3D AMR

Compose the existing cell-local elementary/full-H2/O2 3D reactor with both
levels of the static periodic hierarchy. After each private source phase,
average down the fine state and recover covered coarse temperature. Publish a
complete `R(dt/2)-H(dt)-R(dt/2)` candidate only after both source phases,
hydro, reflux, synchronization, and EOS recovery succeed. Keep the
chemistry-disabled path on the exact qualified hydro dispatch.

Apply the same operation to parent-aligned sparse MPI storage. Permit empty
fine ownership, but require every rank to compare the complete ordered
reaction records, tolerances, integrator, chemistry-enable flag, and existing
hydro metadata before entering or publishing physics. Test serial parity,
rank-dependent valid metadata, invalid source policy, and a failed coarse
source on a zero-fine rank at one, two, four, and eight ranks.

Add full-H2/O2 control/reacting hotspot applications. Check strict CSV schema
and topology, physicality, species closure, reaction activity, five Euler
integrals, H/O/N totals, and exact serial/MPI files. Reject chemistry-enabled
schema-2 checkpoint/restart because its context is incomplete. Defer selected
mechanisms, transport, dynamic topology, physical boundaries, CTU/PPM, 3D EB
AMR, scalable I/O, performance, and external physical validation.

## 0.242.0 selected mechanisms in static 3D AMR

Extract the serial and sparse-MPI application lifecycles from the fixed
frontends, then generate selected frontends that supply complete NASA7,
reaction, transport-provenance, composition, bundle-SHA, and integrator
context. Keep the fixed frontend on the same shared driver and freeze its
full-H2/O2 and characteristic-PLM output hashes.

Before any MPI numerical allocation or file creation, compare the complete
selected context on every rank, including optional presence and bit patterns
of all real fields. Preserve participation by ranks owning no fine planes and
return collective false without mutation for a valid-but-different context.

Qualify a deliberately non-H/O/N two-species generated fixture and generated
full-H2/O2 in serial and at one, two, four, and eight ranks. Require byte-exact
selected serial/MPI files, fixed/selected full-H2/O2 identity, strict generic
output checks, and startup rejection for unsupported transport, checkpoint,
restart, and path aliases. Defer selected persistence, AMR transport, dynamic
topology, physical boundaries, CTU/PPM, 3D EB AMR, scalable I/O, performance,
and external physical validation.

## 0.243.0 selected context-bound static 3D AMR restart

1. Keep fixed schema 2 byte-identical and select schema 3 only when all
   selected context arguments are present.
2. Bind bundle SHA-256, generated integrator, bundle-order composition,
   complete ordered reactions, chemistry tolerances, NASA7 and numerical
   fingerprints, original baseline, history, and both levels.
3. Validate the complete file transactionally through strict EOF before any
   caller state is published; reject incomplete context before target open.
4. Preserve exact communicator-wide selected context consensus and make
   malformed efficiency shapes collective failures before representation
   packing. Include chemistry tolerances in the pre-I/O consensus.
5. Keep the root-formatted checkpoint rank-neutral and demonstrate two-rank to
   one/two/four/eight-rank exact continuation against serial uninterrupted
   output.
6. Freeze selected checkpoint/final hashes, retain fixed hashes, exercise
   cross-schema/context/alias failures, then repeat the full build/test/install
   qualification matrix.

This increment does not add payload authentication, crash-atomic/scalable I/O,
schema migration, AMR molecular transport, dynamic topology, physical coarse
boundaries, 3D EB AMR, performance, or physical/external PeleC validation.

## 0.244.0 molecular transport in static 3D AMR

Reuse the qualified regular 3D mixture transport operator on the static,
strictly interior two-level hierarchy. Select the coarse step from coarse and
fine parabolic limits, advance every fine transport Euler stage with `r^2`
subcycles, and build all fine exterior transport states from coarse start/end
states using time interpolation and limited PLM prolongation. Accumulate
time- and area-averaged fluxes on all six patch faces, then reflux, average
down, and recover temperatures before returning a private candidate.

Protect uncovered coarse species when fine fluxes replace coarse interface
fluxes. Derive one conservative face theta from the coarse state with the old
coarse-face contribution removed, scale the replacement flux, and apply the
equal-and-opposite integrated correction to fine boundary cells. Exercise all
six orientations under an active limiter and require composite conservation.
Compose two hierarchy transport half-steps with chemistry and hydro as a
transactional `R-T-H-T-R` update.

Implement the same operation directly on parent-aligned sparse MPI slabs.
Preserve zero-fine-rank participation, compare gas-transport records and
transport flags before payload communication, reduce deterministic interface
fluxes and limiter values, and retain exact rollback. Qualify serial versus
one-/two-/four-/eight-rank fixed and selected public output, direct adversarial
MPI parity, ratio-three and individual-process unit cases, and invalid theta,
barodiffusion, metadata, and checkpoint/restart paths.

Keep transport checkpoint/restart rejected until a schema binds the transport
database and operator contract. Defer dynamic 3D AMR, boundary-touching fine
patches, physical coarse boundaries, 3D EB AMR, scalable I/O, performance,
detailed-fuel physical validation, and external PeleC field parity.

## 0.245.0 transport checkpoint persistence

Carry the static two-level AMR transport context through an exclusive schema 4
for the fixed public `full_h2o2` case (`chemistry_enabled=.false.`) and schema
5 for selected chemistry (`thermo_model='selected'`,
`chemistry_enabled=.true.`). Both transport-enabled schemas must serialize and
validate the complete ordered gas-transport database, parameter convention,
operator identity, transport controls, cumulative diagnostics, and the
`POST_ACCEPTED_COARSE_STEP` phase after the committed `R-T-H-T-R` update.

Keep schema 2/3 transport-disabled restart byte-compatible and reject all
fixed/selected or transport-enabled/transport-disabled cross-schema fallback.
Read into private level, temperature, clock, baseline, reflux, and transport
diagnostic candidates; require the terminal marker, strict EOF, and close before
publishing any caller state. A failed read must leave all caller targets
unchanged.

Retain root-formatted, rank-neutral I/O with no communicator size or ownership
map. Require all-rank consensus for selected/transport context before payload
communication or root I/O, keep zero-fine-plane ranks in every collective, and
broadcast root status and restored metadata before scattering to current-rank
ownership. Qualify serial and changed-rank MPI continuation exactly for both
fixed and selected transport cases. Formatted replacement remains an
application-level transaction only; crash-atomic replacement, durable commit,
payload authentication, and scalable I/O are not claimed.

The focused 0.245 gates record GNU Debug/Release serial `10/10` each and
GNU/OpenMPI Debug/Release MPI `23/23` each. Known schema-4/schema-5 checkpoint
and output hashes are recorded in
[`validation/0.245.0.md`](validation/0.245.0.md). The full configuration matrix
and clean installed Release qualification are not included in this record.
