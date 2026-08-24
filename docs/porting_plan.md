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
- higher-order StateRedist neighborhood reconstruction and limiting
- thermal, viscous, and catalytic wall boundary conditions

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
