# Completion roadmap

## Definition of complete

PeleF is complete when every PeleC responsibility selected by this project is
either:

1. implemented in independent Fortran with automated numerical evidence; or
2. explicitly excluded with a documented reason and public replacement path.

Completion does not mean transliterating every C++ source line. It means a
runnable reacting-flow solver whose supported dimensions, chemistry,
transport, AMR/EB, parallelism, particles, restart, and output contracts are
precise and reproducible.

## Frozen reference

The comparison baseline is PeleC `development` commit
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`. Its recursive PelePhysics,
AMReX, and SUNDIALS revisions are recorded in
[`references/pelec_baseline.json`](../references/pelec_baseline.json). A future
baseline update must be intentional, must update that manifest, and must state
which numerical signatures changed.

## Workstreams and exit gates

| Order | Workstream | Current state | Exit gate |
|---:|---|---|---|
| 0 | Reproducible project contract | Qualified in `0.201.0` | Version records agree; Debug/Release and MPI configure paths are reproducible; no shipped ELF target requests an executable stack; the upstream manifest is machine checked |
| 1 | Three-dimensional regular and AMR flow | In progress: `0.245.0` adds focused fixed/selected transport checkpoint persistence to the `0.244.0` static periodic two-level 3D AMR PCM/PLM, chemistry, and transport path; `0.217.0` separately exposes stabilized single-level planar 3D EB hydro with optional elementary chemistry, molecular transport, and transactional restart | Full release qualification beyond the focused transport-persistence record; dynamic 3D AMR topology, physical boundaries, and 3D EB AMR are qualified |
| 2 | Production chemistry and transport | Partial: `0.245.0` adds schema-4/schema-5 transport context binding, transactional restart, and focused serial/MPI changed-rank parity to the `0.244.0` `R-T-H-T-R` and `r^2` fine-subcycling path; selected execution and context-bound restart also span the established regular, AMR, EB, and sparse-MPI subsets | Runtime switching, thread-safe and performance qualification of stiff integration, detailed-fuel qualification, and declared transport parity against pinned references remain |
| 3 | Three-dimensional EB | In progress: exact axis-plane geometry, slip-wall fluxes, conservative divergence, PCM Euler, FluxRedist/StateRedist, full-grid hydro/transport limits, and a serial public application with optional elementary chemistry, StateRedist transport, and transactional restart are qualified through `0.217.0` | Cut-cell geometry, wall fluxes, redistribution, reflux, regridding, and restart pass 3D conservation and reference gates |
| 4 | LES models | Not started | Every exposed model has unit, manufactured-flow, and configuration/restart tests |
| 5 | Lagrangian particles and spray | Not started | Injection, interpolation, source coupling, parcel migration, restart, and mass/momentum/energy closure pass serial/MPI gates |
| 6 | Scalable I/O, performance, and accelerators | Not started | Distributed production output/restart, scaling evidence, and the declared CPU/GPU parity envelope pass |
| 7 | Production validation and release | Not started | Representative reacting-flow and spray cases have external validation records; license and release artifacts are complete |

`Partial` means useful verified subsets exist but the responsibility is not
closed. A green structural or rank-parity test is not by itself physical
validation.

## Acceptance rule for each increment

Each increment must contain implementation, a focused unit gate, an
application or parity regression, and documentation of the exact supported
boundary. External comparisons must name the manifest revision, case, initial
state, observable, tolerance, and comparison time. A capability remains open
when only a scaffold, smoke test, or short-time structural case exists.

## Immediate sequence

1. establish a 3D mesh/state/directional-flux core with x/y/z dimensional
   reductions;
2. add a conservative 3D uniform-grid driver and convergence case;
3. extend reactive thermodynamics, chemistry splitting, and molecular
   transport to that path;
4. introduce 3D AMR before 3D EB so reflux and restart invariants are isolated;
5. close mechanism ingestion and production integration before claiming
   detailed-fuel capability;
6. add LES, then particles/spray, on the qualified 3D distributed state.

Milestone `0.202.0` closed the coordinate-core part of step 1. Milestone
`0.203.0` closes step 2 for a serial, single-level, periodic, constant-`gamma`
PCM/SSPRK2 Euler path: x/y/z whole-step reduction is exact, the smooth 3D
entropy wave converges toward first order, and every Euler integral is
conserved to roundoff. Milestone `0.204.0` closes the general-EOS
reacting-state hydrodynamic portion of step 3: elementary/full-H2O2 NASA7
temperature recovery, x/y/z Rusanov/HLLC/PeleC-style dimensional reduction,
and every Euler/species integral now pass. Milestone `0.205.0` closes the
reaction-source portion of step 3 with elementary/full-H2O2 cell reactors,
whole-field transactional Strang composition, elemental conservation, and a
coupled public hotspot gate. Milestone `0.206.0` closes step 3 for the selected
regular-grid transport subset with x/y/z reduction, a full viscous tensor,
thermal/species smoothing, parabolic timestep selection, periodic
conservation, and transactional `R-T-H-T-R` composition. Milestone `0.207.0`
begins step 4 with a serial static two-level, single-patch hydro boundary:
ratio subcycling, time-interpolated coarse boundary states, six-face reflux,
average-down, and composite conservation pass. Milestone `0.208.0` adds a
versioned, physics-fingerprinted hierarchy checkpoint, transactional restart,
and byte-exact continuous/restarted public outputs. Milestone `0.209.0`
assigns coarse/fine x-plane arithmetic to MPI ranks and closes serial versus
1/2/4/8-rank byte parity plus two-to-four/eight-rank restart parity. The
hierarchy is still replicated. Milestone `0.210.0` adds selectable
frozen-composition characteristic PLM to the regular 3D solver and both levels
of the serial static hierarchy, including limited space/time coarse ghost
interpolation, reflux of time-averaged PLM fluxes, and schema-2 restart
fingerprints. Sparse storage, distributed PLM, topology change, physical
coarse boundaries, and CTU/PPM remain open. Milestone `0.211.0` begins the 3D
EB workstream with an exact static, axis-aligned planar geometry on one
Cartesian level. It stores volume and face fractions, cell/face/EB centroids,
EB area, and normal integrals; validates the cellwise discrete geometric
identity; balances general-EOS slip-wall pressure against Cartesian face
pressure; and publishes a transactional PCM Euler candidate. Uniform tangent
flow is exactly invariant and a nonuniform fluid-volume integral is conserved
to roundoff. Oblique/curved or interior face-aligned cuts, small-cell
redistribution, higher-order EB reconstruction, chemistry/transport,
AMR/reflux/regridding, restart, and MPI remain open.

Milestone `0.212.0` adds conservative first-order FluxRedist through the six
open Cartesian neighbor faces. It preserves a uniform active residual,
conserves every fluid-volume-weighted component, and makes a `kappa=0.05`
CFL-0.5 full-grid hydro step physical when the raw cut-cell update produces
negative density.

Milestone `0.213.0` adds a distinct zeroth-order weighted StateRedist route for
the same exact planar geometry. It uses one axis-normal regular receiver,
target-volume and overlap-count partitions, conservative neighborhood
gather/scatter, and transactional NASA7 recovery. Raw, FluxRedist, and
StateRedist now share one hydro residual builder, and both stabilized routes
complete the raw-negative full-grid control. Higher-order StateRedist, a
public EB timestep selector, general geometry, high-order EB reconstruction,
coupled physics, hierarchy operations, restart, MPI ownership, and a public
3D EB application remain open.

Milestone `0.214.0` adds the active-cell full-grid CFL selector and installed
`pelef_reactive_eb_3d` executable. Its validated input requires a stabilized
exact plane, its forward-Euler loop is transactional across repeated steps,
and its CSV exposes geometry plus every reactive conserved field. Independent
StateRedist/FluxRedist public runs and a 480-row checker distinguish physical
wall-normal impulse from the mass, energy, species, and tangential-momentum
invariants. Higher-order integration/reconstruction, configurable chemistry
and transport, general geometry, hierarchy operations, restart, and MPI remain
open.

Milestone `0.215.0` composes optional elementary cell-local chemistry with the
serial planar StateRedist application. The three-dimensional chemistry wrapper
accepts an optional active mask and keeps its candidate state and recovered
temperature private until every plane succeeds. The EB split performs
`R(dt/2)-H(dt)-R(dt/2)`, leaves covered and disabled cells exactly unchanged,
and rejects invalid solver, tolerance, mechanism, geometry, or reaction stages
without publishing a partial result. The application gates mass, total energy,
the two tangential momenta, and H/O/N elemental totals rather than treating
individual species as conserved through reaction. A dedicated unit gate covers
x/y/z masks, rollback, inert hydro parity, species change, and covered identity.
Separate public inert/reacting StateRedist cases compare against the existing
0.214 baseline byte-for-byte and independently check the 480-row geometry,
EOS/closure, active species change, covered identity, fluid invariants, and
element totals. This remains an elementary seven-species, exact-axis-plane,
serial single-level workflow; transport, configurable mechanisms, general
geometry, AMR/reflux/regridding, restart, and MPI chemistry remain open.

Milestone `0.216.0` composes the regular-grid molecular constitutive laws with
the planar 3D EB StateRedist path. Open internal Cartesian faces carry
Newtonian, Fourier, mixture-averaged species, barodiffusive, correction-
velocity, and species-enthalpy fluxes; physical-domain and embedded-wall
transport fluxes are zero. Both transport Euler stages apply the qualified
order-zero StateRedist operation, covered storage is bitwise invariant, and a
physical-inventory limiter bounds outward species transfer. The public step is
the minimum active-state hydro/transport bound and the full transaction is
`R-T-H-T-R`. Direct x/y/z units cover conservation, wall/boundary fluxes,
limiter activation, rollback, and disabled-path parity. Frozen one-step
control, transport, and coupled public cases independently check 480-row
geometry/EOS contracts, effective relative conservation, H/O/N totals,
resolved response, logs, schedules, and full CSV hashes. General isothermal or
no-slip wall transport, FluxRedist transport, higher-order or general geometry,
AMR/reflux/regridding, restart, MPI ownership, and external validation remain
open.

Milestone `0.217.0` adds a serial process boundary to the coupled planar 3D EB
application. Its versioned formatted checkpoint records the full NASA7,
elementary reaction, and gas-transport inputs; immutable mesh, EB, physics, and
numerical settings; state and temperature; initial conserved/L1/element
inventories; and cumulative timestep/transport diagnostics. Reads parse and
validate private candidates through the terminal marker before publishing, so
incompatible, truncated, or EOS-inconsistent files leave every caller field
unchanged. A public chemistry-plus-transport case stops after its first
committed `R-T-H-T-R` step, resumes for at least one further step, and requires
byte-exact final CSV plus diagnostic parity with an uninterrupted run. The
format remains serial, single-level, exact-axis-plane validation I/O; MPI/AMR
checkpoint layouts, general geometry, scalable I/O, and schema conversion are
open.

Milestone `0.218.0` closes the first bounded mechanism-ingestion slice. A
pinned Cantera YAML file and explicit `ohmech` phase are converted at build
time into one deterministic normalized JSON bundle. The importer validates an
ideal-gas phase, NASA7 at 101325 Pa, element balance, Lennard-Jones gas
transport, elementary/three-body/Troe reactions, named colliders, source
ordering, and duplicate flags. The generator emits the thermo, kinetics,
transport, and provenance consumed by the existing full-H2/O2 runtime path.
Live Cantera tests require byte-clean YAML-to-JSON and JSON-to-Fortran
regeneration plus trajectory and multidimensional CFD parity. This does not
close runtime YAML loading, CHEMKIN, NASA9, PLOG/Chebyshev/SRI/Lindemann,
custom reaction orders, mechanisms above 32 species, detailed fuels, or a
production stiff integrator.

Milestone `0.219.0` connects normalized mechanism bundles to the build graph.
An optional CMake selection reads the bundle-declared Fortran interfaces,
regenerates the module in the build tree, compiles it in an isolated library,
and builds an installable probe that executes thermo, reaction, primitive
transport, production-rate, and Jacobian calls. A supplied source YAML is
checked against the bundle basename and SHA-256 during configuration and again
before regeneration. The ordinary production target set remains unchanged;
the test graph compiles a small independent fixture, while a clean
tests-disabled configuration qualifies the pinned full-H2/O2 bundle. The
energy-constrained reactor Jacobian now also passes a centered directional
finite-difference gate. This remains configure-time ABI qualification, not
runtime parsing, fixed-application dispatch, mechanisms above 32 species, new
rate families, detailed-fuel validation, or production stiff integration.

Milestone `0.220.0` turns that selected build artifact into an installed
serial constant-volume application. Exact species-name composition mapping,
finite and tail-complete namelist validation, generator-side element/mass
balance, runtime reaction/mass checks, and common NASA7 range checks precede
output creation. A reactive two-species fixture passes conservation,
activity, and repeat-run byte-identity gates; selecting the pinned full-H2/O2
bundle reproduces the fixed application CSV byte-for-byte. A scheduler gate
also proves an exact final fragment below the adaptive minimum, while the CMake
helper explicitly limits selected builds to complete schema-1 thermo/transport
bundles. After the later `0.223.0` regular-1D boundary, the remaining
production-chemistry exit gate still requires runtime selection and dispatch
through the other CFD/AMR/EB/MPI applications, thread-safe and
performance-qualified stiff integration, detailed fuels, and external
validation.

Milestone `0.221.0` adds the first Phase-7 integration boundary without
changing the native default. An optional dependency-isolated runtime links the
pinned SUNDIALS 7.2.0 official Fortran CVODE, serial-vector, dense-matrix, and
dense-linear-solver modules through a double-precision-only ABI gate. CVODE BDF
advances `N-1` mass fractions, uses the largest initial species for exact
closure, recovers temperature at fixed specific internal energy, and receives
the generalized energy-constrained semi-analytic Jacobian. Invalid trials are
recoverable; failed advances do not publish caller state and require
finalization. Tests cover unavailable-build
rejection before CSV creation, nonfinal dependent-species Jacobian direction,
single-call and cumulative step-limit rollback, singleton rejection,
deterministic two-species output, the short-final-fragment scheduler, full
10-species/29-reaction conservation, static-dependency installation, and
Cantera trajectory/rate parity. This is one serial non-reentrant context, not
runtime mechanism loading, CFD dispatch, sparse chemistry, detailed-fuel
validation, or a performance result.

Milestone `0.222.0` replaces that singleton with an opaque context-first API.
The public handle contains only a private registry slot and generation; a
stable C-interoperable token carrying those values is installed with
`FCVodeSetUserData`, while all noninteroperable mechanism data and solver
resources remain in one of 64 private slots. Calls may be interleaved
sequentially across live contexts. Tests require exact standalone/interleaved
trajectory and statistics identity for states with different closure species,
local failure quarantine while a peer continues, arbitrary and repeated
finalization, stale-copy rejection after slot reuse, and clean capacity
exhaustion. Cleanup reports the first failing vendor destroy routine and then
invalidates the handle; destructor retry is not claimed. This closes serial
multi-context lifecycle qualification only;
threaded/reentrant use, sparse solvers, detailed fuels, CFD dispatch, and
performance/scaling remain open.

Milestone `0.223.0` dispatches a configure-time normalized mechanism through
the existing regular serial 1D reacting-flow algorithm. A separate selected
runtime owns the hydro/configuration/transport modules and links only the
selected generated bundle, avoiding the committed full-H2/O2 module namespace.
The input composition is resolved by exact species name; selected primitive
transport records are validated and converted into the shared mixture model.
A reacting two-species fixture checks activity and invariants, while the pinned
full-H2/O2 selected output is byte-identical to the fixed regular-1D output.
This closes one application boundary only. Runtime mechanism loading, selected
CVODE use inside CFD, 2D/3D/AMR/EB/MPI dispatch, threaded and performance
qualification, detailed fuels, and external validation remain open.

Milestone `0.224.0` extends configure-time dispatch to the periodic regular 3D
application. Common selected-composition and selected-mechanism startup modules
now serve both 1D and 3D front ends. A generic 3D application driver accepts
validated model arrays from the fixed or selected loader, while a private
selected 3D link graph excludes both committed generated mechanisms. The
two-species fixture proves active chemistry and deterministic 4-by-4-by-4
evolution. The selected full-H2/O2 transport hotspot is byte-identical to the
fixed path and retains the independent nonzero transport-response gate. A
short full-H2/O2 characteristic-PLM hotspot additionally proves active selected
chemistry and exact fixed/selected output. This closes periodic single-level
regular 3D only; selected physical boundaries, 2D, AMR, EB, MPI, runtime
loading, CFD CVODE, detailed fuels, and production qualification remain open.

Milestone `0.225.0` dispatches the same configure-time mechanism through the
regular serial 2D application. A shared driver now owns the fixed/selected
simulation, output, invariant, and diagnostic path, while the selected runtime
privately compiles the regular 2D configuration, boundary, CTU, chemistry, and
transport dependencies. The two-species fixture proves active generic
chemistry and deterministic output. Fixed and selected pinned full-H2/O2
executables are byte-identical for uniform chemistry, a nonuniform
characteristic-PLM/CTU chemistry step, and prescribed-species-wall molecular
transport. This closes single-level regular 2D, including the established
physical boundaries. Selected AMR, EB, MPI, runtime loading, CFD CVODE,
detailed fuels, thread safety, and performance qualification remain open.

Milestone `0.226.0` dispatches the same configure-time mechanism through the
serial reactive 1D AMR application. A shared fixed/selected driver owns
two-level, arbitrary-depth, and dynamic-multipatch mode selection, simulation,
composite CSV, conservation, and diagnostics. Optional bundle-order
composition reaches every AMR root initializer without changing fixed callers.
The two-species fixture proves active chemistry, transport, refinement, and
deterministic output. Fixed and selected full-H2/O2 outputs are byte-identical
for active two-level chemistry/transport, three-level characteristic PPM, and
three-patch molecular transport. Selected AMR checkpoint/restart, EB, MPI,
runtime loading, CFD CVODE, detailed fuels, thread safety, and performance
qualification remain open.

Milestone `0.227.0` dispatches the selected bundle through the serial
single-level reactive 2D EB application. A shared fixed/selected driver owns
plane/circle geometry construction, configured boundary assembly, simulation,
deterministic EB CSV, volume-weighted invariants, and diagnostics. Optional
bundle-order composition reaches both regular initialization and reactive
physical-boundary construction, while the generated integrator policy reaches
each active-cell chemistry half-step. The two-species plane fixture proves
active chemistry/transport and deterministic cut-cell output. Fixed and
selected full-H2/O2 outputs are byte-identical for plane chemistry and a
characteristic-PLM circle with isothermal moving no-slip transport. Selected
EB AMR/3D, MPI, restart, runtime loading, CFD CVODE, detailed fuels, thread
safety, and performance qualification remain open.

Milestone `0.228.0` dispatches the selected bundle through the bounded regular
MPI 1D reacting-flow verification application. Fixed and selected front ends
load their own models and call one driver for the uneven 19-cell periodic
decomposition, adaptive `R-T-H-T-R` evolution, collective invariants, ordered
gather, and dynamic-species CSV. The generated integrator policy reaches both
distributed chemistry half-steps and invalid policy returns collectively with
exact state/temperature rollback. The independently ordered two-species
fixture is active and byte-identical across 1/2/4 ranks; selected full-H2/O2
is byte-identical across those ranks, to the fixed MPI executable, and to the
pre-refactor fixed output. This closes neither a general input-driven MPI
application nor selected MPI AMR/EB, restart, scaling, or thread safety.

Milestone `0.229.0` dispatches the selected bundle through the serial static
two-level reactive 2D EB AMR application. Fixed and selected front ends call
one mechanism-independent application driver for configured boundaries,
coarse/fine initialization, fine subcycling, EB reflux, average-down,
deterministic level output, invariants, and diagnostics. The independently
ordered two-species fixture is positive and closed and repeats byte-for-byte
on both levels. The selected pinned full-H2/O2 run changes H2 on both levels
and is byte-identical to the fixed executable under coupled chemistry,
molecular transport, and an isothermal moving no-slip embedded wall. Dynamic,
three-level, multipatch, and restart modes remain rejected until their
selected-mechanism contracts are complete.

Milestone `0.230.0` dispatches the selected bundle through the serial
single-level reactive 3D EB application. Fixed and selected front ends call
one mechanism-independent application driver for exact planar geometry,
adaptive CFL selection, optional `R-T-H-T-R` advancement, redistribution,
diagnostics, and dynamic-species CSV publication. An independently ordered
H2/H fixture is active and repeat-byte identical; a test-only elementary
bundle reproduces the fixed seven-species output byte-for-byte; and the
pinned full-H2/O2 bundle exercises implicit chemistry and all qualified
molecular-transport terms. Selected checkpoint/restart, general geometry and
element diagnostics, EB AMR/MPI 3D, runtime loading, and detailed fuels remain
outside this increment.

Milestone `0.231.0` closes selected checkpoint/stop/restart for that same
serial single-level planar EB 3D boundary. Fixed calls retain byte-identical
schema 1; selected schema 2 records and validates the configure-time bundle
SHA-256, generated integrator policy, and normalized bundle-order composition
before publishing any restored target. A separate selected elementary process
resumes a first-step checkpoint to the exact uninterrupted final CSV and
cumulative diagnostics. The bundle digest is provenance rather than a payload
digest; crash-atomic replacement, tamper evidence, schema migration, general
geometry, and selected AMR/MPI EB 3D persistence remain outside this increment.

Milestone `0.232.0` closes selected checkpoint/stop/restart for the serial
static two-level reactive EB AMR 2D boundary. Fixed calls retain byte-identical
schema 3, while selected schema 4 records and validates the configure-time
bundle SHA-256, generated integrator, and normalized bundle-order composition
before any coarse/fine restart target is published. A pinned full-H2/O2
process resumes its first-step checkpoint to exact uninterrupted coarse and
fine CSV files, and both selected references are byte-identical to the fixed
references. Dynamic, three-level, multipatch, and MPI EB AMR persistence,
payload authentication, crash-atomic/scalable I/O, and physical detailed-fuel
qualification remain outside this increment.

Milestone `0.233.0` closes selected checkpoint/stop/restart for the serial
dynamic exactly-two-level, single-fine-patch reactive EB AMR 2D boundary.
Fixed schema 3 and selected static schema 4 remain byte-identical; selected
dynamic schema 5 binds the same generated context and initial composite
integrals while retaining the actual post-regrid patch and regrid history in
the transactional body. The qualified
full-H2/O2 patch moves from `(2:5,2:5)` to `(5:12,4:12)` before its first-step
checkpoint and resumes to exact uninterrupted coarse and fine CSV files that
also equal the fixed reference, including its cumulative conservation
diagnostic. Selected three-level, multipatch,
dynamic-parent, and MPI EB AMR persistence, payload authentication,
crash-atomic/scalable I/O, and physical detailed-fuel qualification remain
outside this increment.

Milestone `0.234.0` closes selected checkpoint/stop/restart for the serial
static exactly-three-level, one-patch-per-level reactive EB AMR 2D boundary.
Fixed static magic/schema 3 and fixed dynamic magic/schema 4 remain
byte-identical. Selected schema 4 under the static magic binds the generated
context and original composite-integral baseline before the established
three-level body. A pinned full-H2/O2 first-step checkpoint resumes to exact
uninterrupted root, middle, and finest CSV files that also equal the fixed
reference, including the cumulative conservation diagnostic. Selected
dynamic-three-level, multipatch, dynamic-parent, and MPI EB AMR persistence,
payload authentication, crash-atomic/scalable I/O, and physical detailed-fuel
qualification remain outside this increment.

Milestone `0.235.0` closes selected checkpoint/stop/restart for the serial
dynamic exactly-three-level, one-patch-per-level reactive EB AMR 2D boundary,
including dynamic-parent regridding. Fixed dynamic magic/schema 4 remains
byte-identical; selected schema 5 under that magic binds the generated context
and original composite-integral baseline before the established dynamic body.
A pinned full-H2/O2 lifecycle moves both middle and finest patches before its
first-step checkpoint and resumes to exact uninterrupted root, middle, and
finest CSV files that also equal the fixed reference, including the cumulative
conservation diagnostic. Selected multipatch and MPI EB AMR persistence,
payload authentication, crash-atomic/scalable I/O, and physical detailed-fuel
qualification remain outside this increment.

Milestone `0.236.0` closes selected checkpoint/stop/restart for the serial
dynamic two-level multipatch reactive EB AMR 2D boundary. Fixed patch-set
magic/schema 3 remains byte-identical; selected schema 4 under that magic
binds the generated context and original composite-integral baseline before
the established multipatch body. A pinned full-H2/O2 lifecycle commits two
separated children, stops after its first coarse step, and resumes to exact
uninterrupted root and child CSV files that also equal the fixed reference,
including the cumulative conservation diagnostic. Selected MPI AMR/EB
persistence, payload authentication, crash-atomic/scalable I/O, and physical
detailed-fuel qualification remain outside this increment.

Milestone `0.237.0` closes configure-time selected dispatch and rank-neutral
checkpoint/stop/restart for the sparse MPI reactive AMR 1D patch tree. Fixed
schema 1 remains byte-identical; selected schema 2 binds bundle SHA-256,
generated integrator, bundle-order composition, and the original composite
baseline without storing communicator size or ownership. A one-rank
full-H2/O2 checkpoint resumes on two and four ranks to exact uninterrupted
selected and fixed-reference bytes. Selected MPI EB AMR, payload
authentication, crash-atomic/scalable I/O, and physical detailed-fuel
qualification remain outside this increment.

Milestone `0.238.0` closes configure-time selected execution and rank parity
for the sparse MPI reactive EB patch-tree 2D lifecycle. Fixed and selected
front ends share one lifecycle driver; complete selected context is checked
communicator-wide before initialization, bundle-order composition reaches
root and boundary construction, and the generated integrator reaches both
chemistry half-steps on every active level. A four-level dynamic full-H2/O2
case is byte-exact across fixed one rank and selected one, two, and four ranks.
Selected checkpoint/restart remains explicitly unsupported; fixed schema 8 is
frozen against the previous Release binary. Selected sparse MPI EB
persistence, payload authentication, crash-atomic/scalable I/O, and physical
detailed-fuel qualification remain outside this increment.

Milestone `0.239.0` closes selected checkpoint/stop/restart for the sparse MPI
reactive EB patch-tree 2D lifecycle. Fixed schema 8 remains byte-identical;
selected schema 9 binds bundle SHA-256, generated integrator, bundle-order
composition, complete numerical fingerprint, original composite baseline,
clock, counters, regrid history, topology, and all patch fields. A one-rank
full-H2/O2 checkpoint resumes on two and four ranks to exact uninterrupted
fixed and selected output. Payload authentication, crash-atomic/scalable I/O,
runtime loading, fixed-depth MPI modes, and physical detailed-fuel
qualification remain outside this increment.

Milestone `0.240.0` replaces the replicated numerical fields in the public
static MPI 3D AMR application with parent-aligned rank-local coarse and fine x
slabs. Two-plane periodic coarse halos and active-neighbor fine halos support
both PCM and characteristic PLM, including ranks with no fine cells. CFL,
SSPRK2 subcycling, time-interpolated coarse/fine ghosts, six-face reflux,
average-down, and temperature recovery remain local until collective commit.
Ownership, patch, NASA7, timestep, solver, reconstruction, and limiter
metadata are consensus-checked before data-dependent communication. PCM and
PLM outputs are byte-exact against serial at one, two, four, and eight ranks;
two-rank schema-2 checkpoints resume at four and eight ranks exactly. Root
gather/scatter remains rank-neutral compatibility I/O, not scalable I/O.
Dynamic topology, AMR chemistry/transport, physical boundaries, CTU/PPM, 3D
EB AMR, scaling qualification, and external PeleC field parity remain open.

Milestone `0.241.0` composes fixed elementary or full-H2/O2 chemistry with the
static periodic two-level hierarchy. Each source phase advances private
coarse and fine candidates, averages down, and recovers the covered coarse
temperature; the complete step publishes only after
`R(dt/2)-H(dt)-R(dt/2)` succeeds. Sparse ranks with no fine cells participate
in reaction/integrator and final acceptance consensus. The full-H2/O2 hotspot
changes species and temperature, conserves the Euler subset and H/O/N totals,
and is byte-exact between serial and one, two, four, and eight MPI ranks.
Schema-2 chemistry restart is rejected because it lacks reaction and source-
integrator context. Dynamic topology, molecular transport, physical coarse
boundaries, selected mechanisms, CTU/PPM, 3D EB AMR, scalable I/O,
performance, and external field parity remain open.

Milestone `0.242.0` dispatches configure-time selected mechanisms through the
same serial and rank-local sparse-MPI static 3D AMR lifecycle as the fixed
path. Exact-name composition, bundle SHA-256, generated integrator, ordered
NASA7, reaction, and transport provenance are validated before MPI state
allocation. A two-species fixture covers non-H/O/N diagnostics and zero-fine
ranks; generated full-H2/O2 is byte-identical to fixed serial and one-, two-,
four-, and eight-rank output. Unsupported selected transport and persistence
requests fail before output. Context-bound selected restart, AMR molecular
transport, dynamic topology, physical boundaries, CTU/PPM, 3D EB AMR,
scalable I/O, performance, and external field parity remain open.

Milestone `0.243.0` closes selected checkpoint/stop/restart for the serial and
rank-local sparse-MPI static two-level 3D AMR lifecycle. Fixed schema 2 remains
byte-identical; exclusive selected schema 3 binds bundle SHA-256, generated
integrator, bundle-order composition, complete ordered reactions, chemistry
tolerances, NASA7 and numerical fingerprints, original composite baseline,
clock/reflux history, and both levels before transactional publication. The
rank-neutral root-formatted payload contains no communicator size or ownership
map. A serial/two-rank full-H2/O2 first-step checkpoint is byte-identical and
resumes to exact uninterrupted output in serial and on one, two, four, and
eight ranks. The final eight-configuration matrix passes `4366/4366`, and a
tests-disabled clean Release passes its 42/42 installed-ELF audit plus installed
schema-3 and fixed schema-2 changed-rank restart gates. Payload authentication,
crash-atomic/scalable I/O, AMR molecular
transport, dynamic topology, physical boundaries, CTU/PPM, 3D EB AMR,
performance, and external field parity remain open.

Milestone `0.244.0` composes mixture molecular transport with the fixed and
configure-time selected static two-level 3D AMR lifecycle. Each SSPRK2
transport Euler stage advances the fine patch through `r^2` subcycles with
time-interpolated limited-PLM coarse ghosts, accumulates all six time- and
area-averaged boundary fluxes, limits the replacement flux against uncovered
coarse species budgets, applies the conservative fine-side flux correction,
refluxes, and averages down. Serial and parent-aligned sparse MPI use the same
`R-T-H-T-R` transaction; ranks with no fine planes remain valid participants.
Within each build type, the fixed/selected full-H2/O2 public outputs are
byte-identical in serial and at one, two, four, and eight ranks. Transport
persistence remains explicitly rejected. Dynamic topology, physical coarse
boundaries, 3D EB AMR, scalable I/O, performance, detailed-fuel physical
validation, and external field parity remain open. The final matrix passes
`4432/4432` tests; the clean Release install audits `42/42` executables and
repeats installed transport and retained transport-disabled restart parity.

Milestone `0.245.0` persists the static 3D AMR transport context without
changing the `0.244.0` stepping scope. Fixed public `full_h2o2` transport uses
exclusive schema 4 with chemistry disabled; selected chemistry uses exclusive
schema 5 with `thermo_model='selected'`, chemistry enabled, and transport
enabled. Both schemas bind the full ordered gas-transport database, controls,
operator/convention identity, cumulative transport diagnostics, and the
`POST_ACCEPTED_COARSE_STEP` boundary after committed `R-T-H-T-R` stepping.
Reads validate private candidates through strict EOF and close before
publication, while root-formatted serial/MPI I/O remains rank-neutral and
supports changed-rank continuation with all-rank consensus, including ranks
with no fine planes. Schema 2/3 transport-disabled compatibility remains
isolated and cross-schema fallback is rejected.

Focused qualification records GNU Debug/Release serial `10/10` each and
GNU/OpenMPI Debug/Release MPI `23/23` each; schema-4/schema-5 checkpoint and
output SHA-256 values are frozen in
[`validation/0.245.0.md`](validation/0.245.0.md). No `0.245.0` full matrix or
clean installed Release result is reported. Crash-atomic replacement, payload
authentication, scalable I/O, dynamic topology, physical boundaries, 3D EB
AMR, performance, and external PeleC field parity remain open.

The repository currently has no license file. Selecting and adding a license
requires an explicit project-owner decision and is a release blocker, not a
numerical implementation choice.
