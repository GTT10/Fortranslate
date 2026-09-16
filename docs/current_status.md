# Current status: PeleF 0.246.0

This page separates implemented capability, dated reported results, and
independent checks. The preceding integration is
recorded in [PR #102](https://github.com/GTT10/Fortranslate/pull/102). The new
opt-in spray scope and qualification boundary are in
[the spray guide](spray.md) and [0.246.0 validation](validation/0.246.0.md).
Older milestone prose remains historical detail, not a current global status.

## Supported boundary

| Area | Implemented subset | Still open |
| --- | --- | --- |
| Reactive flow and transport | 1D/2D paths; regular 3D NASA7 hydro/chemistry/molecular transport with per-face periodic, prescribed reservoir inflow, outflow, slip-wall and adiabatic/isothermal no-slip boundaries with prescribed tangential wall motion | Characteristic inflow and moving geometry, general 3D boundary/coupling qualification |
| 3D AMR | Two-level hierarchy, with serial/MPI boundary-touching physical patches; serial/MPI temperature-tagged patch movement/resizing with conservative remap, rank-neutral geometry-aware restart and sparse MPI stepping; serial and sparse MPI hydro, reaction and molecular-transport split APIs with physical coarse boundaries and application input/output/restart ([hydro](validation/amr-physical-hydro-20260917.md), [reaction](validation/amr-physical-split-20260917.md), [transport](validation/amr-physical-transport-20260917.md), [application](validation/amr-physical-app-20260917.md)) | Multiple general patches/levels and broader scaling/physical-boundary qualification |
| 3D embedded boundary | Bounded planar, single-level hydro/chemistry/transport/restart paths | General geometry and integrated 3D EB-AMR-MPI qualification |
| Chemistry | Fixed H2/O2 and configure-time selected supported bundles, including Lindemann/Troe/SRI/Tsang falloff, chemically activated rates, PLOG and Chebyshev; official-Fortran CVODE for selected 0D and opt-in regular 3D spray cells; FFCM1_Red methanol, GRI-Mech and Lu n-dodecane fixtures; input/runtime capacity up to 128 species | Diesel ignition/physical validation across conditions and sparse/performance qualification |
| Parallelism/output | Rank-parity and changed-rank restart for qualified subsets; distributed two-level MPI AMR CSV output | Production scaling, broader distributed I/O qualification and other output paths |
| LES and spray | Opt-in regular 3D single-component parcels with per-port liquid identity or shared or per-port multicomponent compositions, density/pressure fits and vapor-state caps, and periodic/reservoir-inflow/outflow/specular parcel boundaries and adiabatic/isothermal tangentially moving gas walls; temperature-dependent liquid density/heat capacity and per-port breakup viscosity/surface tension, selectable Schiller-Naumann or PelePhysics drag with ETAB deformation correction, drag/heating/evaporation with Clausius or bounded Antoine saturation pressure and an optional near-boiling vapor-state cap, and ETAB and KH-RT representative-drop breakup, optional KH shed-mass children with conservative transverse velocities, and bounded computational splitting, conservative two-way coupling and uniform external acceleration with impulse/work ledger, deterministic cone/scheduled injection from point or finite circular solid/hollow outlets with swirl with up to 16 independently configured ports and optional piecewise-linear shared or independent per-port injection-rate histories with optional discharge-based injection speed, and per-port mass- or number-CDF Weibull/normal/lognormal/ChiSquared size distributions, deviatoric Smagorinsky, coupled restart | Physical mixture boiling, general dispersion/collision, particle MPI/AMR/EB, characteristic inlet models, near-wall turbulence, lateral wall-film transport and impact qualification, calibrated properties and experimental validation |

A new regular 3D patch-tree library represents multiple nonoverlapping children
and arbitrary depth, with parent-before-child indexing and a common refinement
ratio at each depth. It reuses existing PCM/restriction kernels and provides
transactional initialization, bottom-up averaging and composite inventories
([design](design_decisions/amr-patch-tree-3d.md),
[evidence](validation/amr-patch-tree-3d-20260917.md)). This is a geometry/state
foundation. Snapshot ghost filling now exchanges peer data across parent boxes,
uses PCM ancestor fallback and follows regular-grid periodic/physical boundary
conventions ([evidence](validation/amr-tree-ghost-3d-20260917.md)). Per-level coarse
endpoint interpolation is also available for ghost donors while retaining exact
current peers ([evidence](validation/amr-tree-temporal-ghost-20260917.md)). PLM
ghost filling, hydro/reflux/subcycle scheduling, chemistry/transport, dynamic
remapping, MPI ownership and application/restart integration remain open.

Sparse MPI AMR has physical ghost layers and two-level PCM/PLM hydro supporting
one coarse column per rank, with serial equality and collective rollback checks
([ghosts](validation/amr-mpi-boundary-halo-20260917.md),
[coarse hydro](validation/amr-mpi-boundary-hydro-20260917.md),
[hierarchy](validation/amr-mpi-physical-hierarchy-20260917.md)). The
[reaction/hydro split](validation/amr-mpi-physical-reaction-20260917.md) also forwards
physical boundary data and budgets. Physical molecular transport and the full
reaction/transport/hydro split are also connected
([evidence](validation/amr-mpi-physical-transport-20260917.md)). MPI application
inputs, wall/inlet timestep limits and rank-neutral boundary-ledger restart are
also connected ([evidence](validation/amr-mpi-physical-app-20260917.md)). Fixed elementary
and full-H2/O2 reacting AMR checkpoint/restart now retains canonical reaction,
composition and integrator context for serial/MPI applications
([evidence](validation/amr-fixed-chemistry-restart-20260917.md)).

Boundary-touching AMR now has a physical flux-correction API that replaces
covered domain-surface fluxes in the boundary ledger and corrects only uncovered
interior coarse cells ([evidence](validation/amr-boundary-reflux-20260917.md)).
Serial PCM/PLM hydro now connects fine physical ghosts, subcycling, boundary
ledger replacement and average-down for such patches
([evidence](validation/amr-boundary-touching-hydro-20260917.md)). Cell-local reaction,
reaction/hydro splitting and the state-based CFL helper now accept these patches
([evidence](validation/amr-boundary-touching-reaction-20260917.md)). Molecular
transport and the full reaction/transport/hydro split are connected, including
fine wall heat/work fluxes ([evidence](validation/amr-boundary-touching-transport-20260917.md)).
Serial static application inputs, fine wall timestep bounds and checkpoint/restart
are connected ([evidence](validation/amr-boundary-touching-app-20260917.md)).
Serial temperature-driven regridding now supports boundary-touching patches with
explicit physical tagging and geometry-aware restart
([evidence](validation/amr-boundary-regrid-20260917.md)). Sparse MPI distribution,
scatter/gather and cell-local reaction accept boundary-touching geometry
([evidence](validation/amr-mpi-boundary-distribution-20260917.md)). MPI PCM/PLM
hydro and reaction/hydro splitting now connect fine physical ghosts and covered
surface flux replacement ([evidence](validation/amr-mpi-boundary-touching-hydro-20260917.md)).
MPI molecular transport and the full reaction/transport/hydro split now also
connect fine wall heat/work and physical flux replacement
([evidence](validation/amr-mpi-boundary-touching-transport-20260917.md)).
MPI static application inputs now include fine wall timestep limits and
rank-neutral restart for boundary-touching geometry
([evidence](validation/amr-mpi-boundary-touching-app-20260917.md)).
MPI physical-boundary temperature tagging, conservative remap and changed-rank
restart are also connected
([evidence](validation/amr-mpi-boundary-regrid-20260917.md)). Temperature tagging
now scans local slabs and reduces only tag bounds; unchanged patches skip full
hierarchy gathering ([evidence](validation/amr-mpi-distributed-tagging-20260917.md)).
Conservative remapping now runs within fixed coarse slabs using neighboring
coarse halos and local fine overlap reuse
([evidence](validation/amr-mpi-distributed-remap-20260917.md)). It no longer
gathers full fields during regridding. An explicit repartition API now transfers
coarse/fine states and temperatures directly between ranks, and remapping
preserves the chosen coarse partition
([evidence](validation/amr-mpi-repartition-20260917.md)). Optional cell-count
balancing is connected after initialization/restart and geometry changes
([evidence](validation/amr-mpi-load-balance-20260917.md)). General patches/levels,
measured performance/scaling and broader physical qualification remain open. MPI AMR
CSV output now uses direct collective file writes, and final inventory, synchronization
and composition diagnostics use distributed reductions
([evidence](validation/amr-mpi-distributed-csv-20260917.md)). Analytic initial
conditions are also generated on local coarse/fine slabs, with local average-down
and distributed initial inventories
([evidence](validation/amr-mpi-local-initialization-20260917.md)). Checkpoint writes
now serialize the existing format directly from local slabs after collective
metadata, temperature and synchronization checks
([evidence](validation/amr-mpi-checkpoint-output-20260917.md)). Restart now reads
canonical checkpoint cell rows directly into local candidate slabs, then commits
geometry, fields and diagnostics only after collective checks
([evidence](validation/amr-mpi-checkpoint-input-20260917.md)). Metadata remains
replicated; reformatted cell records use a bounded-memory streaming fallback.
Production I/O scaling is not yet measured.

Ahamed wall-impact classification, splash velocities and azimuthal mass
fractions, four-child construction and deposited cylindrical film patches are
connected to the regular coupled-step API with wall budgets, stationary film
heating/evaporation with per-material constant or cubic thermal conductivity,
CLI inputs/output and coupled checkpoint persistence
[evidence](validation/wall-film-conductivity-20260916.md). Lateral film transport
and physical impact qualification remain open. The gas-side film skin now uses
the pinned half-temperature/half-composition rule
([correction evidence](validation/wall-film-skin-20260916.md)).

KH shed-mass accumulation and four-child release now connect to regular CFD
and restart ([evidence](validation/khrt-children-20260916.md)); computational
splitting now partitions the KH inventories
([evidence](validation/khrt-splitting-20260916.md)). Stationary wall impact
now consumes crossed parcels before KH release, with coupled film/lineage restart
([evidence](validation/khrt-wall-20260916.md)); wet-impact correlation limits remain.

The mixed-liquid API now covers mixture properties, component transport,
regular-grid feedback, component-specific boundary budgets, representative
breakup/splitting and composition checkpoint persistence
([evidence](validation/spray-mixture-grid-20260916.md)). The combined
hydro/transport/LES/chemistry driver now retains composition and rolls it back
with all other state ([evidence](validation/spray-mixture-reactive-20260916.md)).
Mixture CLI inputs, shared/independent scheduled injection, component CSVs and
restart compatibility are connected
([evidence](validation/spray-mixture-cli-20260916.md)). Trace-species depletion and
minimum-diameter removal now return their complete inventories to gas
([evidence](validation/spray-mixture-terminal-20260916.md)). Mixed KH children now inherit composition with conservative transverse velocities
and retain their state through splitting and restart
([evidence](validation/spray-mixture-kh-20260916.md)). The mixture impact event API now partitions component masses into splash and film
with wall budgets ([evidence](validation/spray-mixture-wall-event-20260916.md));
its CLI connection is now available. Planar crossing now uses mixture properties
and the reference mass-weighted impact boiling temperature
([evidence](validation/spray-mixture-wall-crossing-20260916.md)). A supplied film-exchange API now
updates composition, temperature and height with gas/wall energy closure
([evidence](validation/spray-mixture-film-exchange-20260916.md)). Mixed film
heat/evaporation rates and a frozen-rate update are now available
([evidence](validation/spray-mixture-film-rates-20260916.md)). The checkpoint API
now retains film component masses with transactional validation
([evidence](validation/spray-mixture-film-checkpoint-20260916.md)); application wiring is connected. A transactional mixed-film grid feedback API now
conserves species and gas/film/wall energy and removes terminal patches
([evidence](validation/spray-mixture-film-grid-20260916.md)). The mixed particle
coupler now connects wall crossing, splash/deposition and remaining-interval film
evolution in one transaction ([evidence](validation/spray-mixture-wall-coupled-20260916.md)).
The reactive driver now retains mixed film state across both spray half-steps and
rolls it back after late failure ([evidence](validation/spray-mixture-wall-reactive-20260916.md)).
CLI wall inputs, film component output, conservation totals and restart are connected
([evidence](validation/spray-mixture-wall-cli-20260916.md)). Active trace-film
components now transfer their full residual inventories to gas
([evidence](validation/spray-mixture-film-trace-20260916.md)). Per-port compositions now drive initial and
scheduled injection with component budgets and restart identity
([evidence](validation/spray-port-mixture-20260916.md)). Component density/heat-capacity fits, Antoine
pressure and Clausius surface caps are available in the CLI and restart context
([evidence](validation/spray-mixture-properties-20260916.md)).

Conservation, rollback, byte-parity, and short-time H2/O2 tests do not establish
physical ignition/spray validation or complete PeleC equivalence. Choosing a
license remains an owner decision.

Particle MPI groundwork now includes Cartesian cell ownership and stable send-buffer
partitioning with component/KH state preservation
([evidence](validation/spray-partition-20260916.md)). MPI migration and receive-side ID ordering now pass 2/4-rank checks
([evidence](validation/spray-mpi-migration-20260916.md)). Distributed gas/film
feedback and application integration remain open. Collective 64-bit ID reservation
now rejects overflow and inconsistent rank histories
([evidence](validation/spray-mpi-ids-20260916.md)). Collective computational splitting
now uses these IDs while preserving component and KH inventories
([evidence](validation/spray-mpi-split-20260916.md)). Collective cone injection now
reserves IDs and retains mixed component budgets with rollback
([evidence](validation/spray-mpi-injection-20260916.md)). Distributed field sampling
now requests only interpolation cells from owner ranks and matches the serial
stencil ([evidence](validation/spray-mpi-sampling-20260916.md)). Owner-routed
feedback deposition now matches serial stencil accumulation
([evidence](validation/spray-mpi-deposition-20260916.md)). Gas feedback now recovers
owner-block temperatures and commits conserved state and temperature collectively
([evidence](validation/spray-mpi-gas-20260916.md)). Periodic mixed-particle drag/heating/evaporation now uses owner sampling, a
common substep, collective feedback retry and end-of-interval migration with
whole-interval rollback ([evidence](validation/spray-mpi-coupling-20260916.md)).
Representative ETAB and KH-RT breakup now run inside the MPI mixed-particle
candidate update with component-inventory preservation
([evidence](validation/spray-mpi-breakup-20260917.md)). Breakup CLI settings and
deformation/RT checkpoint state now support changed-rank continuation
([evidence](validation/spray-mpi-breakup-app-20260917.md)). A collective KH shedding
API now constructs children with unique IDs, inherited composition and a gas-energy
ledger ([evidence](validation/spray-mpi-kh-release-20260917.md)). MPI coupled
feedback, input and shedding/lineage restart are now connected, including ID
rollback on rejected updates
([evidence](validation/spray-mpi-kh-coupled-20260917.md)). Periodic x-slab SGS now
shares the serial Smagorinsky face kernel, with MPI SSPRK2, coupled half-steps,
CLI coefficients and restart identity
([evidence](validation/spray-mpi-les-20260917.md)). Stationary mixed wall-film
patches can now migrate to Cartesian owners while retaining geometry and component
inventories ([evidence](validation/spray-mpi-film-migration-20260917.md)). An MPI
film-coupling API now connects heating/evaporation to owner-routed gas feedback
with a global wall-energy ledger and whole-interval rollback
([evidence](validation/spray-mpi-film-coupled-20260917.md)). MPI checkpoints now
collect and redistribute film/component state and the wall budget across rank
changes ([evidence](validation/spray-mpi-film-restart-20260917.md)). A collective
mixed-wall impact event API now reserves splash IDs and commits parent, child,
film and wall-budget updates transactionally
([evidence](validation/spray-mpi-impact-20260917.md)). Distributed film-footprint
queries now sum overlapping film heights across owners
([evidence](validation/spray-mpi-film-query-20260917.md)). An MPI planar-crossing
wrapper now computes mixture impact conditions, queries film thickness and resolves
wall events ([evidence](validation/spray-mpi-wall-crossing-20260917.md)); integration
with nonperiodic gas/particle stepping and the application remains open.
A two-layer MPI ghost API now combines physical wall/inlet/outlet conditions with
neighbor exchange, including transverse corners
([evidence](validation/spray-mpi-boundary-ghost-20260917.md)). MPI hydro now uses
these boundaries in both SSPRK2 stages and returns a global physical-face flux
ledger ([evidence](validation/spray-mpi-boundary-hydro-20260917.md)).
MPI molecular transport now shares the serial physical-boundary flux kernel,
including isothermal conduction and tangential wall work, with a global boundary
ledger ([evidence](validation/spray-mpi-boundary-transport-20260917.md)).
MPI SGS now uses the same physical-boundary ghost and wall-heat rules and
returns a global surface ledger ([evidence](validation/spray-mpi-boundary-les-20260917.md)).
The MPI particle-coupling API now supports open-face removal, stationary wall
reflection and nonperiodic owner migration with a replicated boundary ledger
([evidence](validation/spray-mpi-particle-boundary-20260917.md)). Mixed wall-impact/film
creation now participates in the coupled driver transaction
([evidence](validation/spray-mpi-impact-coupled-20260917.md)). Application impact/film
inputs, per-rank film outputs, total inventories and splash-lineage restart now
connect these models ([evidence](validation/spray-mpi-film-app-20260917.md));
physical qualification remains open.
Physical boundary data now pass through the complete MPI particle/chemistry/
transport/SGS/hydro split, with a summed boundary ledger and whole-step rollback
([evidence](validation/spray-mpi-boundary-split-20260917.md)). Application inputs now
expose physical faces and wall data; history and rank-neutral restart retain the
cumulative boundary budget ([evidence](validation/spray-mpi-boundary-app-20260917.md)).
Existing distributed films can now heat/evaporate in symmetric half-steps around
the complete coupled update, with film/gas rollback and wall-budget accumulation
([evidence](validation/spray-mpi-film-split-20260917.md)).
Periodic particle/reaction splitting now connects cell-local chemistry between
particle half-steps with collective whole-interval rollback
([evidence](validation/spray-mpi-reaction-20260916.md)). Distributed wall/breakup
events remain open. Periodic
x-slab gas advection now exchanges two halo layers and uses characteristic PLM /
SSPRK2, including optional integration between reaction half-steps
([evidence](validation/spray-mpi-hydro-20260917.md)). Molecular transport now
exchanges gas and species-limiter halos and connects as half-steps around hydro
([evidence](validation/spray-mpi-transport-20260917.md)). The fixed-mechanism
`pelef_mpi_spray_3d` executable now connects initial two-component cone injection,
global time-step selection/retry, coupled stepping, budget history and rank-local
CSV output ([evidence](validation/spray-mpi-application-20260917.md)). Scheduled
independent multiport injection now merges up to 16 port schedules with deterministic
simultaneous-event order and changed-rank restart
([evidence](validation/spray-mpi-multiport-20260917.md)). Independent per-port rate
tables now override the shared profile with complete restart identity
([evidence](validation/spray-mpi-port-profiles-20260917.md)). Finite-radius nozzle,
solid/hollow/swirl and size-distribution controls now reach the MPI injector
([evidence](validation/spray-mpi-nozzle-20260917.md)). Shared piecewise-linear injection
profiles now provide per-event integrated mass, including zero-dose event handling
and restart identity ([evidence](validation/spray-mpi-profile-20260917.md)). A shared MPI application
now also serves `pelef_mpi_spray_3d_selected`, with 21-species FFCM1_Red execution
and restart evidence ([validation](validation/spray-selected-mpi-20260917.md)).
The selected case also exposed and corrected roundoff-negative final-species mass
closure ([correction](validation/reactive-zero-species-closure-20260917.md)). Repeated timed cone events now clip time steps, retain injected budgets and
resume from the saved event count ([evidence](validation/spray-mpi-schedule-20260917.md)).
MPI spray CLI save/restart
now preserves model identity and budgets, with focused same-rank and 2↔4-rank
continuation evidence ([validation](validation/spray-mpi-restart-application-20260917.md)).
The checkpoint API now collects mixed particles and owner gas blocks into a root
file and redistributes them with a different communicator size
([evidence](validation/spray-mpi-checkpoint-20260917.md)).
The checkpoint API now preserves a 64-bit ID high-water mark with strict validation,
ready for MPI reservation state ([evidence](validation/spray-id-checkpoint-20260917.md)).

The spray CVODE branch now builds after removal of misplaced wall-film setup;
focused native/CVODE chemistry tests pass
([repair evidence](validation/spray-cvode-build-repair-20260916.md)).

## Repository consolidation implemented

The README is an entry point rather than a chronological implementation log.
Contributor instructions, AGENTS.md, the documentation index and PR template
separate current capability, historical evidence and future work.

The 2613-line root CMake and 9138-line test registration were split into five
root target modules and nineteen test modules, preserving include order and
directory scope. Generated target graphs, install rules and CTest commands
matched in seven clean configurations; see [the layout record](validation/build-layout-20260909.md)
and [CMake navigation](../cmake/README.md).

Hosted workflows cover source provenance/regeneration, complete serial and MPI
suites, selected 0D CVODE with pinned SUNDIALS 7.2.0 and Cantera 3.2.0, and a
separate tests-disabled Release installation. Source identity, dependency
versions, JUnit and logs are retained. An implemented workflow is not itself a
passed qualification; completed results must be associated with its tested SHA.

Historical draft #11 was reviewed and closed as superseded: its 17 numerical,
mechanism, application/case and verification paths were byte-identical at the
later integrated full-chemistry commit. Its unique decision note is preserved
[in the historical record](history/pr11-implicit-chemistry-decision.md).
No historical branch was deleted. The audited ancestor-only bulk cleanup was
blocked by the connection's safety restriction, not performed silently.

## Evidence ledger

The saved checkpoint is `509d07929d52d9e570a6af1f6ec7e5836bd4fe5d`.
Initial PR CI failed because its vendored YAML description lacked two spaces.
Restoring the exact Cantera package bytes repaired source provenance without
changing the mechanism bundle or generated Fortran; see [the integrity repair](validation/checkpoint-integrity-20260908.md).

A later targeted numerical review found a false-success return in NASA7 molar
conversion overflow guards. Commit `bb2980e52c2e6a4bfed353dd280fe1f0e21f4893`
corrects the failure flag without changing physical formulas or tolerances.
New cp/enthalpy/entropy rejection cases fail on the parent implementation and
pass after the fix; seven focused GNU14 Debug tests passed. See [red/green evidence](validation/nasa7-conversion-status-20260909.md).
This correction is intentionally separate from the behavior-preserving layout.

**Reported historical local qualification:** the [2026-09-08 handoff](handoff-pelef-0.245.0-20260908.md)
reports 4570/4570 across Debug 524, Release 524, Debug+Cantera 529, MPI Debug 694,
MPI Release 694, CVODE Debug 533, CVODE Debug+Cantera 539, and CVODE Release 533.
It also reports tests-disabled 42/42 ELF, separate tests-enabled 23/23 installed
restart checks, and installed application/checker/byte comparisons of 32/32,
8/8 and 12/12. These belong to that frozen worktree and its local evidence;
the raw server directories were not retrieved here. Interrupted 227-test and
excluded literal-prefix exploratory results are not successful qualification.

**Independent source-integrity qualification:** commit
`566010ffc77a6de25d452d084ff61fa2bb533d95` passed
[run 34208154353](https://github.com/GTT10/Fortranslate/actions/runs/34208154353),
including 14 checker regressions, exact hashing and full bundle regeneration.
Final serial/MPI/CVODE/clean-install results belong to PR #102 and their own
source identities. They must not be replaced by older source successes.

The earlier focused 10/10 serial and 23/23 MPI results remain in
[validation/0.245.0.md](validation/0.245.0.md). Historical evidence is preserved,
not overwritten or added together as though all configurations were identical.

See [regular 3D boundary inputs](../cases/reactive_boundaries_3d/README.md) for the new
hydro/chemistry boundary scope and flux accounting.

## Remaining work

[The roadmap](completion_roadmap.md) distinguishes further solver development
from repository consolidation. General 3D AMR/EB, production detailed-fuel stiff CFD,
full LES/spray qualification and external physical validation remain genuine engineering work,
not features implicitly completed by housekeeping or large test counts.
