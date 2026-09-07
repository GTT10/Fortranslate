# Fortranslate / PeleF

PeleF is an independent Modern Fortran reimplementation of selected numerical algorithms and capabilities from PeleC. It is not a mechanical C++ translation and is not an official Pele Suite project.

Reference implementation: `Pele-Suite/PeleC:development` at
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`; recursive revisions are frozen in
[`references/pelec_baseline.json`](references/pelec_baseline.json).

## Current capability

The `0.245.0` milestone contains the released `0.244.0` capability baseline
and adds focused static-3D-AMR transport checkpoint persistence. The released
`0.244.0` baseline contains seventeen default serial
executables, ten
serial configure-time selected-mechanism executables, eleven fixed optional
MPI executables, four configure-time selected-mechanism MPI executables, and
runnable serial and sparse-MPI one-dimensional
reactive AMR applications with solution-driven dynamic regridding and
molecular transport. It now includes a conservative uniform-grid `pelef3d`
application with periodic x/y/z flux divergence, multidimensional CFL
selection, transactional SSPRK2 stepping, and a smooth entropy-wave
convergence and conservation gate. The new `pelef_reactive_3d` path applies
the same conservative boundary to an elementary or full-H2O2 NASA7 mixture,
recovers temperature transactionally at every SSPRK2 stage, and conserves all
Euler and species components. It now optionally composes cell-local chemistry
as `R(dt/2)-H(dt)-R(dt/2)` over one private 3D candidate, with exact rollback
of state and temperature on any reaction or hydro rejection. The public
periodic hotspot exercises both chemistry and hydrodynamic response while an
independent checker enforces Euler and H/O/N elemental conservation. The same
3D application now optionally adds periodic Newtonian viscosity, Fourier
conduction, mixture-averaged species diffusion, barodiffusion, correction
velocity, and species-enthalpy flux. It selects the minimum hydrodynamic and
parabolic stable steps and publishes the complete `R-T-H-T-R` update only
after every stage succeeds. A full-H2O2 transport/control pair checks the
public configuration, transport database, conservation, physicality, and
resolved thermal/species response. The regular 3D hydro path now also offers
frozen-composition characteristic PLM in x, y, and z with MC or minmod
limiting. Its diagonal smooth-wave gate observes orders `2.02` and `2.46`.
The new `pelef_amr_reactive_3d` path adds
a static, strictly interior, ratio-two fine patch with coarse-step SSPRK2
fluxes, two fine substeps, time-interpolated coarse boundary states, six-face
reflux, and conservative average-down. Its full-H2O2 public entropy wave
checks composite conservation and coarse/fine synchronization. With transport
disabled, the hierarchy writes and transactionally reads a versioned,
self-describing checkpoint
containing both levels, recovered temperatures, time/step history,
conservation baseline, reflux history, mesh/solver settings, and the complete
NASA7 species fingerprint. A stop/restart public case is byte-identical to an
uninterrupted run at both levels. The MPI 3D AMR application now stores only
each rank's coarse x slab and the fine x planes whose parents it owns. Two
periodic coarse ghost planes and active-rank fine halos support both PCM and
characteristic PLM without materializing either global level during a step.
Ranks with no fine cells remain valid participants. Reflux, average-down,
temperature recovery, and candidate publication remain owner-local and
collectively transactional. Its PCM and PLM 1/2/4/8-rank outputs are
byte-identical to the serial results. A two-rank root checkpoint can be
resumed on four or eight ranks with the same exact final fields; root gather
is retained only as rank-neutral compatibility I/O. Fine-patch faces use two
coarse ghost layers filled by time interpolation and limited spatial
prolongation; reflux still consumes the SSPRK2 time-averaged high-order face
fluxes. On the public entropy wave, fine-grid density L1 decreases from
`1.2481e-3` for PCM to `4.3313e-4` for PLM, while composite conservation
remains below `1e-15`.
The static two-level 3D AMR path now also supports fixed elementary or
full-H2/O2 cell-local chemistry in serial and rank-local sparse MPI runs.
Each step is a hierarchy-wide transactional `R(dt/2)-H(dt)-R(dt/2)` update;
both chemistry phases average down the fine state and recover covered coarse
temperatures before publication. Serial and one-, two-, four-, and eight-rank
hotspot outputs are byte-identical, including ranks that own no fine cells.
The application gates the five Euler integrals and H/O/N totals separately
from reacting species activity. Chemistry-enabled schema-2 checkpoint/restart
is rejected because that schema does not persist reaction or integrator
context.
The static two-level 3D AMR path now also dispatches configure-time selected
mechanisms in serial and rank-local sparse MPI execution. A shared lifecycle
driver accepts the generated NASA7, ordered reaction, transport-provenance,
composition, bundle-SHA, and integrator context. Before root state
initialization, MPI ranks compare that complete context bit-for-bit. The
two-species fixture remains exact at one, two, and four ranks, including ranks
with no fine cells; the full-H2/O2 selected path is byte-identical to both the
fixed executable and one-/two-/four-/eight-rank MPI output. Selected
checkpoint/restart now uses an exclusive schema 3 that binds the bundle SHA,
generated integrator, bundle-order composition, complete ordered reactions,
chemistry tolerances, NASA7 records, numerical fingerprint, baseline, and both
levels before transactional publication. The root-formatted file contains no
rank count or ownership map, so a two-rank full-H2/O2 checkpoint resumes
exactly on one, two, four, and eight ranks. The same fixed and selected static
hierarchy now accepts mixture molecular transport in serial and sparse MPI
execution. Each transport half-step is SSPRK2, advances the fine patch with
`r^2` subcycles per Euler stage, accumulates all six time- and area-averaged
interface fluxes, refluxes, and averages down before publication. A common
coarse/fine interface limiter protects uncovered coarse-cell species budgets
and applies the equal-and-opposite correction to fine boundary cells. The
full hierarchy composes chemistry, transport, hydro, transport, and chemistry
as one transaction. Fixed/selected serial and one-/two-/four-/eight-rank MPI
transport outputs are byte-identical within each build type; transport
checkpoint/restart was explicitly rejected in the released `0.244.0` path
because neither checkpoint schema recorded its operator contract.

The current `0.245.0` increment adds exclusive schema 4 for the fixed public
`full_h2o2` case (`chemistry_enabled=.false.`, `transport_enabled=.true.`) and
schema 5 for selected chemistry (`thermo_model='selected'`,
`chemistry_enabled=.true.`, `transport_enabled=.true.`). Both schemas bind the
complete ordered gas-transport database, transport controls, cumulative
transport diagnostics, operator identity, and `POST_ACCEPTED_COARSE_STEP`
phase. Reads validate private candidates through strict end-of-stream and
publish only transactionally; root-formatted serial/MPI checkpoints remain
rank-neutral, so exact changed-rank continuation is covered while the old
schema 2/3 transport-disabled compatibility contracts remain isolated. Focused
0.245 validation records serial GNU Debug/Release `10/10` each and GNU/OpenMPI Debug/
Release MPI `23/23` each; no full matrix or clean-install result is claimed.
See [`docs/validation/0.245.0.md`](docs/validation/0.245.0.md).

The final `0.244.0` matrix passes `4432/4432` tests across eight build
configurations. An isolated tests-disabled Release install audits `42/42`
executables and repeats fixed/selected serial/MPI transport and retained
transport-disabled restart parity using only installed programs.
The first 3D embedded-boundary increment now constructs exact x-, y-, or
z-normal planar cut-cell metrics, validates the discrete aperture/normal
identity, evaluates general-EOS slip-wall pressure fluxes, and advances a
transactional PCM Euler candidate with zero-gradient outer faces. Uniform
tangential flow is invariant under all three Riemann solvers and a nonuniform
case conserves every fluid-volume-weighted component to `2.78e-17`. A new
six-face conservative FluxRedist stage now blends every cut-cell right-hand
side with its open-face neighborhood and returns the excess with exact
fluid-volume weighting. For a volume fraction `0.05`, a CFL-0.5 full-grid
step that makes the raw cut density negative instead completes with density
`0.8725725`. The same planar kernel now also provides zeroth-order weighted
StateRedist on the provisional conserved state. Its normal-directed
neighborhood reaches target volume `0.5`, preserves every
fluid-volume-weighted component, and completes that control step with cut
density `0.2078409091`. Both stabilized routes are transactional through
NASA7 temperature recovery. The installed `pelef_reactive_eb_3d` application
now selects either route, advances repeated full-grid CFL steps, and writes
cell type, volume fraction, fluid centroid, primitive state, and all species
fields to deterministic CSV. Its public two-step density-sheet cases conserve
mass, energy, species, and tangential momentum to `2.13e-15`; wall-normal
momentum is reported separately as physical wall exchange. Arbitrary or curved
geometry, face-aligned interior cuts, higher-order StateRedist, AMR,
and MPI EB remain open. The `0.215.0` increment optionally composes
the existing three-dimensional cell-local elementary chemistry with this
StateRedist workflow as a transactional `R(dt/2)-H(dt)-R(dt/2)` update.
Chemistry is masked to active cells, preserves covered storage bit-for-bit, and
recovers each active temperature from the conserved state. The paired inert
output is byte-identical to the frozen 0.214 StateRedist output, while the
reacting run changes active species and preserves mass, total energy,
tangential momentum, and H/O/N totals to the reported roundoff gates. The
`0.216.0` increment adds StateRedist-only planar molecular transport with
Newtonian viscosity, Fourier conduction, mixture-averaged species diffusion,
barodiffusion, correction velocity, and species enthalpy flux. It selects the
minimum hydro/transport step and publishes the coupled `R-T-H-T-R` candidate
only after every stage succeeds. Outer-domain and embedded-wall transport
fluxes are zero, covered storage remains bitwise fixed, and independent
control/transport/coupled outputs freeze roundoff conservation plus resolved
transport and chemistry responses. The `0.217.0` increment adds a versioned,
transactional checkpoint for this coupled planar path. It fingerprints the
complete thermodynamic, reaction, transport, geometry, and numerical contract,
preserves cumulative diagnostics, and requires byte-exact final CSV parity
after a one-step stop/restart. The `0.218.0` increment adds deterministic,
build-time ingestion of the pinned Cantera `h2o2.yaml` `ohmech` phase. The
normalized bundle preserves NASA7 thermodynamics, ideal-gas Lennard-Jones
transport, element compositions, all 29 source reaction indices, duplicate
flags, named colliders, third-body efficiencies, and the Troe falloff record.
The generated Fortran now supplies the existing full-H2/O2 thermo, kinetics,
and transport runtime path, while CTest regenerates both JSON and Fortran and
compares them byte-for-byte. The `0.219.0` increment adds an optional CMake
selection boundary for any normalized bundle within that schema. It generates
the declared Fortran module in the build tree, compiles it as a separate
library, and builds an installable probe that executes its thermo, reaction,
transport, production-rate, and Jacobian interfaces. An optional source YAML
is checked against the bundle's SHA-256 and basename at configure and build
time. This is build-time qualification, not runtime CFD dispatch. General wall
or FluxRedist transport, general
geometry, AMR, MPI chemistry/transport, scalable I/O, runtime YAML loading,
CHEMKIN, and arbitrary detailed mechanisms remain open.
The `0.220.0` increment added the installed `pelef0d_selected` application to
that configure-time boundary. It maps initial mole fractions by exact bundle
species name, loads the generated NASA7 and reaction tables, rejects invalid
or mass-unbalanced inputs before creating output, and advances the existing
adaptive implicit constant-volume reactor. A two-species fixture is
conservative and byte-deterministic, and the selected pinned H2/O2 trajectory
is byte-identical to the fixed full-H2/O2 application. The selected build
requires a complete schema-1 NASA7/transport bundle; legacy kinetics-only JSON
remains generator-only. The adaptive minimum permits only an unreduced shorter
fragment needed to land on a scheduled final time. CSV `wdot_*` values are in
`kmol/(m^3 s)`. This is a serial 0D application for normalized bundles of
2--32 species, not runtime mechanism loading, CFD application dispatch, or
physical validation of an arbitrary mechanism.
The `0.221.0` increment adds an optional SUNDIALS `7.2.0` CVODE BDF backend to
that same selected application through the official Fortran interface. The
native implicit backend remains the default and is SUNDIALS-independent.
Configuration requires the official static targets and a double-precision
SUNDIALS ABI; single- and extended-precision installations are rejected.
CVODE advances a private reduced mass-fraction state, reconstructs the largest
initial species as the closure species, recovers temperature from fixed
specific internal energy, and uses the energy-constrained semi-analytic
Jacobian with a dense linear solver. Failed advances do not modify caller
state; that initial adapter permitted one serial, non-reentrant context per
process. The two-species fixture is byte-deterministic, the short-final-step
contract is retained, and the pinned 10-species/29-reaction H2/O2 trajectory
passes conservation and Cantera 3.2.0 gates. This does not qualify selected
chemistry in CFD, mechanisms above 32 species, detailed fuels, concurrent
reactors, sparse solvers, or performance.
The `0.222.0` increment replaces that singleton adapter with an opaque,
generation-checked handle API. Up to 64 CVODE contexts may be live in one
process and advanced in any sequentially interleaved order. A stable
C-interoperable `(slot, generation)` token is passed through CVODE user data;
the noninteroperable Fortran mechanism state and all SUNDIALS resources remain
in a private registry. Failure quarantine, cumulative step budgets,
statistics, and cleanup are context-local. Finalizing one context does not
disturb another, copied stale handles cannot release a reused slot, and
capacity exhaustion leaves the requested handle empty. The interleaved
two-state gate is exactly identical to the corresponding standalone
trajectories and solver statistics. This is serial multi-context support, not
thread-safe or reentrant execution, performance qualification, CFD chemistry
dispatch, or broader mechanism validation.
The `0.223.0` increment brings the configure-time selected mechanism into an
installed regular-grid `pelef_reactive_1d_selected` application. It maps an
input composition by exact species name, constructs the mixture-transport
database from the selected bundle, and advances the existing conservative
reactive hydro, chemistry, and molecular-transport path without linking the
committed fixed-H2/O2 generated module. A reactive two-species fixture checks
activity, conservation, deterministic output, and invalid-name rejection.
For the pinned ten-species/29-reaction H2/O2 bundle, the selected and fixed
regular-1D CSV files are byte-identical. This is configure-time dispatch for
one serial regular 1D application, not runtime mechanism loading, CVODE-backed
CFD chemistry, 2D/3D/AMR/EB/MPI dispatch, thread-safety or performance
qualification, detailed-fuel validation, or external PeleC field parity.
The `0.224.0` increment extends that isolated mechanism boundary to the
installed periodic `pelef_reactive_3d_selected` application. A shared generic
3D application driver now receives already validated thermo, reaction,
transport, and exact-name composition arrays from either the fixed or selected
front end. A two-species 4-by-4-by-4 reactor proves active chemistry and
repeat-run byte identity. A short full-H2/O2 characteristic-PLM hotspot proves
active selected chemistry and is byte-identical to the fixed front end. The
selected pinned full-H2/O2 transport hotspot separately exercises viscosity,
conduction, mixture diffusion, barodiffusion, correction velocity, species
enthalpy, and SSPRK2 hydro and is also byte-identical to the fixed 3D output.
This remains serial, periodic, regular-grid, configure-time dispatch; selected
2D/AMR/EB/MPI, runtime loading, CFD CVODE, detailed-fuel validation, thread
safety, and performance remain open.
The `0.225.0` increment adds the installed
`pelef_reactive_2d_selected` application without copying the regular 2D time
loop. Fixed and selected front ends supply validated model arrays and a
bundle-order composition to one generic 2D application driver. The isolated
selected 2D runtime supports the existing periodic and physical-boundary PCM,
characteristic-PLM/PPM, CTU, native chemistry, and molecular-transport paths
for 2--32 species. An independent two-species fixture proves active chemistry
and repeat-byte identity. Matching full-H2/O2 uniform chemistry,
characteristic-PLM hotspot chemistry, and prescribed-species-wall transport
cases are byte-identical between fixed and selected executables. This closes
serial single-level regular 2D configure-time dispatch, including its existing
physical boundaries; selected AMR/EB/MPI, runtime loading, CFD CVODE,
detailed-fuel validation, thread safety, and performance remain open.
The `0.226.0` increment adds the installed
`pelef_amr_reactive_1d_selected` application. Fixed and selected front ends
now call one mechanism-independent AMR driver, while an isolated selected AMR
runtime compiles the existing two-level, arbitrary-depth, and dynamic
multipatch 1D hierarchy modules without either committed generated mechanism.
The independently ordered two-species fixture proves active chemistry,
molecular transport, regridding, and repeat-byte identity. Matching pinned
full-H2/O2 cases are byte-identical for active two-level chemistry/transport,
three-level characteristic PPM, and three-patch molecular transport. This is
also the first selected-flow release to require each complete schema-1 bundle
to declare `chemistry_integrator = "explicit"` or `"implicit"`. The generated
constant is threaded through selected regular 1D/2D/3D and AMR 1D chemistry,
so an arbitrary ten-species bundle is no longer classified by species count.
Fixed front ends retain their established policy. A boundary-touching outflow
case additionally exercises three levels, characteristic PPM, and hybrid
WENO7-Z. This is
serial 1D configure-time selected AMR dispatch without checkpoint/restart;
selected EB/MPI, runtime loading, CFD CVODE, detailed-fuel validation, thread
safety, and performance remain open.
The `0.227.0` increment adds the installed
`pelef_reactive_eb_2d_selected` application. Fixed and selected front ends
load their own models and call one mechanism-independent single-level EB
application driver. The isolated selected runtime composes the existing plane
or circle geometry, PCM or characteristic-PLM hydro, native explicit or
implicit chemistry, molecular transport, StateRedist, and embedded-wall
controls for schema-1 bundles containing 2--32 species. An independently
ordered two-species plane fixture proves active chemistry, transport, cut-cell
masking, and repeat-byte identity. Pinned full-H2/O2 plane chemistry and
isothermal moving no-slip circle-transport cases are byte-identical between
fixed and selected executables. This closes only configure-time selected
serial single-level 2D EB dispatch; selected EB AMR/3D, MPI, restart, runtime
loading, CFD CVODE, detailed-fuel validation, thread safety, and performance
remain open. All four outer-domain boundaries remain outflow-only in this EB
application.
The `0.228.0` increment adds the installed
`pelef_mpi_reactive_1d_selected` verification application. Fixed and selected
front ends load their own models and call one mechanism-independent MPI
driver over the established 19-cell uneven periodic decomposition. A private
selected MPI runtime contains the generic domain, hydro, native chemistry,
and molecular-transport modules without linking `pelef_core` or either
committed mechanism. The generated explicit/implicit policy reaches both
distributed Strang half-steps. The two-species fixture is reactive and
rank-invariant for 1/2/4 ranks; the pinned full-H2/O2 selected output is
byte-identical across ranks and to the fixed MPI output. This is a bounded
configure-time regular-1D MPI verification boundary, not selected MPI AMR/EB,
general input-driven MPI flow, restart, runtime loading, CFD CVODE,
detailed-fuel validation, thread safety, or performance qualification.
The `0.229.0` increment adds the installed
`pelef_reactive_eb_amr_2d_selected` application for the established static
two-level serial EB hierarchy. Fixed and selected front ends share geometry,
boundary construction, subcycled `R-T-H-T-R` advancement, reflux,
average-down, deterministic coarse/fine CSV, and diagnostics. An independently
ordered two-species fixture proves geometry, positivity, closure, and exact
repeat output on both levels. The pinned full-H2/O2 selected run uses the
bundle-declared implicit integrator, changes H2 on both levels, and is
byte-identical to the fixed executable with chemistry, molecular transport,
and an isothermal moving no-slip embedded wall. Dynamic, three-level,
multipatch, checkpoint/restart, MPI EB AMR, runtime loading, CFD CVODE,
detailed-fuel validation, thread safety, and performance remain open.
The `0.230.0` increment adds the installed
`pelef_reactive_eb_3d_selected` application for the established serial,
single-level, exact planar 3D EB path. Fixed and selected front ends share
geometry, CFL selection, `R-T-H-T-R` advancement, StateRedist, diagnostics,
and dynamic-species CSV output. A reversed-order two-species fixture proves
active chemistry and repeat-byte identity. A generated seven-species bundle
using the fixed elementary kinetics, thermodynamics, and transport produces
byte-identical fixed/selected Debug and Release CSV files, while the pinned
full-H2/O2 bundle exercises implicit chemistry and all transport processes.
Chemistry-enabled species outside the qualified H/O/N/Ar set remain rejected
before output.
The `0.231.0` increment qualifies selected checkpoint/stop/restart for that
same serial planar EB 3D path. Fixed calls retain byte-identical schema-1
checkpoints. Selected calls use schema 2 and persist the configure-time bundle
SHA-256, generated explicit/implicit integrator policy, and normalized
bundle-order composition before the existing complete model, geometry, state,
and diagnostic records. Fixed and selected readers reject the other schema;
selected bundle, policy, or composition mismatch fails transactionally. A
separate-process selected elementary case stops after step one and resumes to
the exact uninterrupted 1152-row CSV. AMR/MPI EB 3D, general geometry,
runtime loading, arbitrary-element diagnostics, payload tamper evidence,
detailed-fuel validation, thread safety, and performance remain open.
The `0.232.0` increment qualifies checkpoint/stop/restart for the selected
serial static two-level reactive EB AMR 2D path. Fixed calls retain their
byte-identical schema-3 stream. Selected calls use exclusive schema 4 and
place the configure-time bundle SHA-256, generated integrator policy, species
count, and normalized bundle-order composition ahead of the unchanged
two-level configuration, patch, coarse/fine state, and terminal records.
Cross-schema reads and changed SHA, policy, or composition fail before any
restart target is published. A full-H2/O2 selected process stops after coarse
step one and resumes to the exact uninterrupted coarse and fine CSV files;
the uninterrupted selected files are also byte-identical to the fixed path.
Checkpoint/restart remains limited to static exactly-two-level selected EB
AMR: dynamic, three-level, multipatch, MPI EB AMR, runtime loading, scalable
I/O, payload tamper evidence, detailed-fuel validation, thread safety, and
performance remain open.
The `0.233.0` increment extends that selected persistence boundary to serial,
dynamic, exactly-two-level reactive EB AMR 2D with one fine patch. Fixed calls
retain byte-identical schema 3, static selected calls retain byte-identical
schema 4, and dynamic selected calls use exclusive schema 5. Schema 5 binds
the same configure-time bundle context and stores the initial composite
integrals needed to retain cumulative conservation diagnostics, while the
established body records the dynamic policy, actual refined patch, regrid
count, clock, and coarse/fine fields. The qualified full-H2/O2 hierarchy moves
its fine patch
from coarse bounds `(2:5,2:5)` to `(5:12,4:12)` before the first-step
checkpoint, then resumes in a separate process to byte-exact uninterrupted
coarse and fine CSV files that also match the fixed reference. Static/dynamic
cross-schema reads fail transactionally, and the uninterrupted/restarted
maximum composite conservation error is exact. Three-level, multipatch, dynamic
At that milestone, three-level, multipatch, dynamic-parent regridding, MPI EB
AMR, runtime loading, scalable or crash-atomic I/O, payload authentication,
detailed-fuel validation, thread safety, and performance remained open.
The `0.234.0` increment extends selected persistence to the established serial,
static, exactly-three-level reactive EB AMR 2D hierarchy with one patch per
level. Fixed static three-level checkpoints remain byte-identical schema 3;
selected static calls use schema 4 under the same checkpoint magic and add the
complete configure-time bundle context plus the original composite-integral
baseline. Restart validates all three private state/temperature candidates and
the terminal marker before publishing any caller target. A pinned full-H2/O2
process stops after coarse step one at `6.2673103087327447e-9 s` and resumes to
the exact uninterrupted root, middle, and finest CSV bytes at `1.0e-8 s`;
those selected reference files also equal the fixed reference at every level,
and the cumulative conservation diagnostic is textually identical. Selected
dynamic three-level, multipatch, dynamic-parent, MPI EB AMR, runtime loading,
scalable or crash-atomic I/O, payload authentication, detailed-fuel
validation, thread safety, and performance remain open.
The `0.235.0` increment extends selected persistence to the established serial,
dynamic, exactly-three-level reactive EB AMR 2D hierarchy with one patch per
level, including the dynamic-parent lifecycle. Fixed dynamic three-level
checkpoints remain byte-identical schema 4; selected dynamic calls use schema 5
under the same dynamic magic and prepend the complete configure-time bundle
context plus the original composite-integral baseline to the established body.
The qualified full-H2/O2 process moves the middle patch from `(2:11,2:11)` to
`(5:10,3:10)` and the finest patch from `(6:9,6:9)` to `(3:10,3:14)` before
its first-step checkpoint, then resumes in another process to the exact
uninterrupted root, middle, and finest CSV bytes at `1.0e-8 s`. Those selected
files also equal the fixed reference at every level, and the cumulative
conservation diagnostic is textually identical. Selected multipatch and MPI
EB AMR persistence, runtime loading, scalable or crash-atomic I/O, payload
authentication, detailed-fuel validation, thread safety, and performance
remain open.
The `0.236.0` increment extends selected persistence to the established serial,
dynamic two-level multipatch reactive EB AMR 2D lifecycle. The fixed patch-set
checkpoint remains byte-identical schema 3; selected calls use schema 4 under
the same magic and prepend the configure-time bundle context plus the original
composite-integral baseline. The qualified full-H2/O2 process commits two
separated 10-by-10 child patches over a 14-by-14 root, stops after coarse step
one at `4.5125253966820509e-9 s`, and resumes in another process to step three
at `1.0e-8 s`. Fixed/selected reference and uninterrupted/restarted root and
both child CSV files are byte-exact, including the cumulative conservation
diagnostic. Selected MPI AMR/EB persistence, runtime loading, scalable or
crash-atomic I/O, payload authentication, detailed-fuel validation, thread
safety, and performance remain open.
The `0.237.0` increment adds the installed
`pelef_mpi_amr_reactive_1d_selected` application and context-bound,
rank-neutral restart to the existing sparse MPI AMR 1D patch tree. Fixed
schema 1 remains byte-identical; selected schema 2 records the generated
bundle SHA-256, integrator, bundle-order composition, and original composite
baseline while excluding communicator size and owner assignments. A pinned
full-H2/O2 process stops on one rank and resumes on two and four ranks to the
exact uninterrupted selected and fixed-reference CSV bytes. Bundle,
integrator, composition, species-order, baseline, rank-consensus, and path-
alias failures publish no output. Density, raw and near-floor species,
closure, temperature, and geometry corruptions are also rejected before
restart publication. Selected MPI EB AMR, runtime loading,
scalable or crash-atomic I/O, payload authentication, detailed-fuel
validation, thread safety, and performance remain open.
The `0.238.0` increment adds the installed
`pelef_mpi_reactive_eb_patch_tree_2d_selected` application on the existing
owner-local sparse MPI reactive EB patch-tree lifecycle. Fixed and selected
front ends now share one mechanism-independent driver. The generated bundle-
order composition initializes the root and boundary states, while the
generated explicit/implicit policy reaches both chemistry half-steps on every
active level. A four-level dynamically regridded full-H2/O2 case produces
byte-identical fixed one-rank and selected one-, two-, and four-rank composite
CSV output. Complete selected context is rank-consensus checked before field
initialization. Selected checkpoint/restart remains intentionally rejected,
including checkpoint schedules, while the fixed schema-8 Release checkpoint
and stopped CSV remain byte-identical to `0.237.0`. An independent two-species
fixture target also executes the generated explicit policy and validates its
four-level sparse EB output.
The `0.239.0` increment admits selected checkpoint/stop/restart on that shared
sparse MPI EB lifecycle through exclusive schema 9 while preserving fixed
schema 8 byte-for-byte. Schema 9 binds the generated bundle SHA-256,
integrator, bundle-order composition, complete numerical fingerprint, original
composite baseline, clock, operator counters, regrid history, rank-neutral
topology, and all patch fields. A one-rank full-H2/O2 checkpoint resumes on
two and four ranks to byte-identical uninterrupted fixed and selected output.
Context, schema, fingerprint, geometry, baseline, density, raw or near-floor
species, closure, temperature, terminal, truncation, trailing-content, and all
six input/output/checkpoint/restart alias failures publish no output. The
MPI adapter also rejects rank-dependent selected metadata presence
collectively before root I/O or gather. The formatted root I/O is not
crash-atomic, scalable, or payload-authenticated.
The preceding
`0.202.0` coordinate core
also covers x/y/z directional flux consistency for ideal-gas,
passive-multispecies, and reacting states. The `0.201.0` work makes the
completion roadmap and upstream
reference machine-readable, rejects an incompatible `mpi_f08` toolchain during
configuration, keeps build/runtime/documented versions synchronized, and
requires non-executable stacks for shipped ELF targets. The sparse MPI driver
can write an intermediate
patch-tree checkpoint and restart it with a different MPI rank count. The
public sparse-MPI arbitrary-depth 2D EB application now qualifies that same
process boundary: a two-rank run writes a four-level checkpoint using one
ownership weighting, then four- and eight-rank processes restart it using
another weighting and reproduce an uninterrupted one-rank reference. The
fresh sparse-MPI 2D EB path now constructs its numerical root state on the
single owning rank and moves those arrays directly into sparse storage;
non-owners never allocate a root state or temperature field. The public serial
and sparse-MPI patch-tree checkpoints now store a schema-8
physics, mesh, EB, and regridding fingerprint and reject incompatible restart
inputs while still permitting changed final time, output/checkpoint schedule,
MPI rank count, and ownership weighting. Multilevel EB conservation closure
now distributes residuals only to active, unrefined parent cells in a local
three-by-three support band around each direct coarse/fine interface; it no
longer perturbs every unrefined cell of the parent patch. Cut-cell molecular
transport now adds first-order normal Fourier heat transfer for isothermal
embedded walls and Newtonian momentum transfer plus wall work for no-slip
embedded walls. The validated embedded-wall control travels with the existing
boundary set, so the same kernel is used by single-level, fixed-depth,
multipatch, arbitrary-depth, serial, and sparse-MPI transport paths; the
default remains adiabatic, free slip, and species impermeable. The
single-level public EB application now reads `embedded_wall_kind`,
`embedded_wall_thermal`, `embedded_wall_temperature`, and the three-component
`embedded_wall_velocity` from `&embedded_boundary`. It rejects isothermal or
no-slip selections unless the matching thermal or viscous transport operator
is enabled. Fixed-depth and arbitrary-depth AMR applications now apply the
same configured wall. Their checkpoint contracts record the wall kind,
thermal mode, temperature, velocity, transport switches, and transport CFL,
and reject incompatible restarts transactionally. The EB AMR library also
exposes conservative MC-limited linear prolongation. Regular parents use
Cartesian child offsets. Cut parents form slopes between fluid-volume
centroids, remove the volume-fraction-weighted mean child offset, and limit
every conserved component to the active coarse-neighbor envelope. Covered or
topology-mismatched parents retain PCM, and an inadmissible linear candidate
retries transactionally with PCM. The
fixed-depth public AMR configuration now selects `pcm` or `linear` through
`prolongation_method` in `&eb_amr`; the selection is used by two-level,
dynamic-regrid, separated sibling-patch, and three-level initialization. The
public hot-wall transport case selects `linear`. Because the selection is not
yet stored in checkpoint identity, fixed-depth checkpoint/restart and the
serial or sparse-MPI arbitrary-depth patch tree explicitly remain PCM-only at
the `0.184.0` boundary. In `0.185.0`, every fixed-depth checkpoint advances to
schema 3 and the shared arbitrary-depth fingerprint advances to schema 4 to
record the method. Linear initialization, dynamic regridding, checkpointing,
and restart are therefore active in the public two-level, multipatch,
three-level, serial arbitrary-depth, and sparse-MPI arbitrary-depth cases. In
`0.186.0`, the same configured linear path therefore preserves each cut
parent's volume-weighted conserved average while retaining nonconstant
active-child states next to the embedded boundary. In `0.187.0`, cut-parent
slopes use a connected multidimensional least-squares fit of the active coarse
fluid centroids, including a minimum-norm rank-one fit where the EB leaves only
one resolved direction. In `0.188.0`, a rank-deficient 3-by-3 fit grows to a
face-connected 5-by-5 stencil before accepting a rank-one fallback, recovering
smooth two-dimensional variation around narrow or turning fluid paths without
crossing covered geometry. In `0.189.0`, the fixed three-level library can
transactionally relocate or resize its root-to-middle patch from root
temperature tags. It first restricts the old finest state into the middle,
regrids the middle while retaining same-resolution overlap, rebuilds a valid
interior finest patch, and publishes all three levels only after a composite
conservation check. In `0.190.0`, `dynamic_parent_regridding` connects that
transaction to the public three-level initialization and periodic schedule.
Dynamic checkpoint schema 4 stores both actual refined patch descriptors and
rejects a restart whose parent-regridding policy differs. In `0.191.0`, the
arbitrary-depth 2D EB tree qualifies recursively tagged outflow-boundary
children in both public applications. Every populated level reaches the same
physical side, physical-side reflux is omitted, and sparse 1/2/4/8-rank output
retains serial field identity. In `0.192.0`, that boundary-touching topology
also survives a selected-root checkpoint, independent restart, and a 2-to-4/8
rank ownership change with uninterrupted field parity. In `0.193.0`, the same
fresh and restarted boundary trees execute thermal conduction through the
recursive `R-T-H-T-R` schedule in serial and at every qualified sparse-MPI
rank count. In `0.194.0`, those lifecycles also enable Newtonian viscosity,
mixture-averaged species diffusion, barodiffusion, correction velocity, and
species enthalpy flux. In `0.195.0`, both chemistry half-steps are active too,
so fresh and restarted boundary trees execute the complete transactional
`R-T-H-T-R` schedule. In `0.196.0`, schema-5 patch-tree checkpoints also retain
the cumulative minimum transport limiter across serial and sparse-MPI restart,
instead of resetting that diagnostic at the continuation boundary. In
`0.197.0`, schema-6 checkpoints retain the original composite-integral
baseline too, so the reported final conservation error continues to cover the
complete logical run instead of only its post-restart suffix. In `0.198.0`,
schema-7 checkpoints also retain cumulative chemistry, transport, and hydro
patch advances for every configured AMR level. The fixed-capacity vectors
preserve dormant deeper-level history when regridding temporarily shrinks the
tree, and sparse MPI restores the same global counters after a rank-count or
ownership change. The checkpoint also retains the cumulative number of
scheduled regrid evaluations and tagged
cells in `0.199.0`. Schema-8 restart therefore reports the same AMR adaptation
history as an uninterrupted run, including after sparse ownership and rank
count changes. In `0.200.0`, the public full-physics case carries two separated
temperature features so the runtime tree branches while one branch continues
to touch the x-upper physical boundary. Fresh serial/sparse execution,
independent serial restart, and two-rank checkpoint continuation at four and
eight ranks now retain both branches with identity-keyed field parity. The tree
can also write one composite CSV containing
every leaf cell exactly once; sparse MPI gathers numerical nodes only to a
selected writer root and reports completion collectively. A dedicated serial
`pelef_reactive_eb_patch_tree_2d` application now reads the established 2D
reactive/EB/AMR namelists, builds up to a configured depth from temperature
tags, advances the complete `R-T-H-T-R` physics clock, regrids periodically,
and writes that composite CSV. The
serial two-dimensional EB AMR driver can create, move, resize, remove, and
re-create one fine rectangle from temperature-gradient tags while preserving
its composite conserved state, and can compose active-cell chemistry with the
two-level or root-only EB hydrodynamic path. Formatted checkpoints preserve
that lifecycle state for a later serial restart. With multipatch mode enabled,
the same application clusters disconnected tags into multiple fine rectangles,
selects a set-wide CFL limit, advances their reactive hydrodynamics and
chemistry transactionally with one flux register per child, periodically
rebuilds the set, and writes one CSV per active child. A dedicated formatted
checkpoint preserves the complete ordered patch set for transactional serial
restart. A configured static fine rectangle may also meet an outflow physical
boundary: its physical-side exterior state is extrapolated from the current
fine boundary cell while the remaining coarse/fine sides retain coarse-time
interpolation and conservative reflux. Temperature-gradient tagging and both
single- and multipatch planners can now create such outflow-side rectangles
dynamically, including one-sided gradient detection on root boundary cells.
The single-patch, separated sibling-patch, and strictly nested three-level
lifecycles also accept mixture molecular transport. Every parent/child pair
advances SSPRK2 transport with ratio subcycling, time-interpolated parent
exterior states, an independent open-area diffusive flux register, reflux,
and hierarchy-wide average-down inside the transactional `R-T-H-T-R`
composition.
The MPI EB AMR bridge executes chemistry, owner-tiled finite-band root hydro
and transport, and fine ratio subcycles on deterministic physics owners. Fine
transport owners assemble patch-plus-two start, uncorrected-end, and current
corrected state/temperature directly from intersecting root tile owners, along
with child-intersecting x/y flux fragments. They reconstruct the four-edge
exterior context, assemble the compact coarse register, perform reflux locally,
and return corrected support directly to those tile owners. Fine state remains
on its owner, and no complete root state or flux array enters the child route.
Collective input consensus and deferred publication preserve exact all-rank
rollback.
One outer MPI transaction now composes those owner operators as
`R-T-H-T-R`. It publishes the root, children, limiter minimum, and chemistry,
hydro, and transport counters only after every stage succeeds; a rejection in
the middle hydro stage discards valid reaction and transport prefixes.
The 2D EB bridge also has a rank-local sparse payload container. Each rank
allocates conserved state and temperature only for its owned root tiles and
children; an explicit materialization boundary reconstructs the replicated
hierarchy when a legacy operator or output path still requires it. A second
materialization boundary gathers each root tile and child only to a
caller-selected root rank. Non-root ranks keep the complete output unallocated,
and one packed point-to-point message is sent per remote entity. This is the
field-sparse foundation for checkpoint and output adapters. Sparse MPI wrappers
now connect that boundary to the existing formatted multipatch checkpoint and
root/child CSV writers. Only the selected writer allocates the complete fields
or touches the files, and its I/O result is returned collectively. The inverse
restart boundary reads the formatted checkpoint only on a selected root and
sends each root tile or child directly to its current owner. Non-root complete
read buffers stay unallocated and no numerical-field broadcast is used. The
restart caller now supplies a replicated geometry-only child descriptor; it no
longer carries replicated child state or temperature fields. Compatibility
wrappers retain the former full patch-set API for existing callers.
Separately, a geometry-only 2D EB patch-tree topology now represents any number
of refinement relations. Each level may branch across multiple ordered parents
and children, and a complete replacement plan commits only after parent links,
EB geometry, patch nesting, and sibling separation validate.
Chemistry now runs directly on those sparse owner allocations. Root tiles and
children are reacted locally with covered cells masked; only post-reaction
average-down communicates numerical state. Each child owner sends one
coarse-footprint restriction buffer only to the distinct root tile owners that
intersect it; unrelated ranks receive nothing. Those owners recover temperature
and apply the restriction locally. Chemistry no longer materializes a complete
hierarchy or broadcasts child restrictions.
The complete sparse `R-T-H-T-R` transaction now remains on sparse owners from
input through commit. Both chemistry half-steps, both SSPRK2 transport
half-steps, and hydro call their direct sparse entrypoints; no complete child
hierarchy is materialized between operators. Counts and the limiter minimum are
still published only with the final outer commit.
Hydro also has a direct sparse entrypoint. Each root tile owner receives only
the neighboring row fragments required for a six-row halo, advances its own
bounded EB band, and retains its stage-start, stage-end, current corrected
state/temperature, and uniquely owned flux rows locally. Each child owner
assembles patch-plus-two start/end/corrected state directly from intersecting
tile owners and extracts the four-edge interpolation context locally. The same
tile owners route their x/y interface-flux fragments directly to the child.
After compact register accumulation, ratio subcycling, and reflux, corrected
support returns directly to those tile owners in deterministic child order.
The final corrected root rows commit locally. No complete hydro state,
temperature, or flux result is assembled on a root physics owner, and there is
no post-compute tile-result or final-scatter traffic.
SSPRK2 molecular transport now follows the same sparse ownership boundary.
For each Euler stage, a root tile owner receives only the neighboring row
fragments needed for its six-row transport/StateRedist guard, advances that
target EB band, and retains its owned stage-start, stage-end, corrected state,
temperature, and unique flux rows locally. A periodic y-boundary tile uses a
boundary-anchored cyclic band built from two contiguous source-row fragments.
One extra row isolates
the required six-row guard from the deliberate internal gap; a small root that
cannot hold that guard uses the complete root band as its local compute band.
No post-compute complete root state, temperature, or flux bundle is assembled.
Each fine child assembles patch-plus-two stage-start,
uncorrected-end, and current corrected state/temperature directly from
intersecting root tile owners. Those owners also send their coarse x/y flux
fragments. The child extracts its exterior context, assembles its register,
performs ratio subcycling, fine-flux accumulation, and reflux locally, then
returns corrected fragments directly to the same root tile owners. The final
SSPRK2 root blend and EOS
temperature recovery are cell-local on those owners, so that step needs no root
gather or scatter. For cut-interface closure, tile owners sum their physical-
boundary flux contributions and combine one conserved `nvar` vector with
`MPI_Allreduce`; unrefined-cell corrections remain tile-local. A late child
failure leaves every sparse allocation bitwise unchanged and publishes zero
work or transfers.
EB flux registers now store only the patch-plus-one-cell correction support.
For sparse transport, each root tile owner retains its computed x-flux rows and
uniquely owned y-faces, then sends only child-intersecting fragments directly
to that child owner. The child assembles patch-local face rectangles and
accumulates its coarse register locally; the root physics owner no longer
builds or sends that register. The established complete-root accumulation call
is a compatibility wrapper over the same compact kernel. Remote child owners
allocate no complete root state, temperature, or flux field; fine state and
coarse correction no longer make a root-physics-owner round trip for reflux.
Exterior-context extraction also accepts globally indexed patch-plus-one
coarse start/end state and temperature support. The complete-root entrypoint
is a wrapper over that support kernel. Sparse MPI uses it on the child owner
after direct tile-fragment assembly without changing interpolation or EOS
recovery.
The sparse hierarchy also selects its own stable coarse interval. Every root
tile evaluates its EB hydro and molecular-transport limits directly on its
owner, while fine children do the same and scale their stable fine steps by the
refinement ratio. One communicator-minimum reduction selects the interval with
zero root-field traffic. Invalid owner state rejects collectively with zero
published dt and transfer count.
The same sparse hierarchy now owns a public full-physics time loop. It
recomputes that distributed stable interval before every `R-T-H-T-R` step,
clips the last interval exactly to the requested target time, and publishes
time, step, operator-count, limiter, and timestep-traffic diagnostics only
after each complete step commits. A later step limit or physics rejection
preserves the already committed prefix and its matching diagnostics.
An explicit sparse regrid can now change the ordered child-patch topology,
reapply the serial EB overlap/prolongation rules, rebuild deterministic
subcycle-weighted ownership, and return to one-copy sparse storage. The whole
  distribution/template/state replacement is transactional. Regrid now remains
  field-sparse: old child owners restrict directly to intersecting root owners,
  root tiles initialize each distinct new child owner, and retained overlap
  moves directly from its old owner to its new owner.
The sparse public clock can also evaluate temperature-gradient tags at a
caller-selected accepted-step cadence. Only the root physics owner constructs
the patch plan; compact topology metadata is broadcast, while the caller's
geometry callback rebuilds each planned EB child. Physics and any due regrid
commit as one transaction, so a geometry or topology failure leaves the step,
clock, hierarchy, counters, and published traffic unchanged.
The EB transfer foundation also accepts one strictly nested three-level
hierarchy, computes its composite integral without double counting, and
average-downs its generic or reactive state from the deepest level to the root
as one rollback-safe transaction.
The same hierarchy can advance reactive EB hydrodynamics recursively with one
root update, ratio-subcycled middle and finest updates, an independent flux
register at each interface, and final deepest-first synchronization. When the
finest interface crosses the embedded boundary, a multilevel conservation
closure returns the measured mass, total-energy, and species residual to
uncovered active middle cells before the outer reflux.
The three-level driver now composes active-cell reaction half-steps on the
root, middle, and finest meshes around that recursive hydro transaction. A
second deepest-first reactive average-down makes the post-chemistry hierarchy
authoritative before all three levels are published together.
With `three_level_enabled`, the public EB AMR application constructs that
static hierarchy from two nested namelist rectangles, selects a root timestep
from all three active-cell CFL limits, advances to the requested final time,
and writes separate root, middle, and finest CSV fields. A dedicated
three-level checkpoint preserves all three conserved and temperature fields,
accepted time and step accounting, and the full nested topology for a
transactional serial restart. With dynamic regridding enabled, the middle
level remains fixed while EB-aware temperature-gradient tags can move and
resize the finest patch transactionally after initialization and accepted
steps. A distinct dynamic three-level checkpoint stores that committed finest
rectangle, regrid count and compatibility controls so a split run resumes the
actual hierarchy and cadence rather than the configured seed.
The single-level reactive EB application now also composes mixture viscosity,
thermal conduction, species diffusion, and barodiffusion with chemistry and
hydrodynamics. Diffusive Cartesian-face fluxes are interpolated to EB face
centroids, weighted by open area, stabilized with EB-aware species limiting
and StateRedist, and use an adiabatic slip/impermeable embedded-wall closure.

### `pelef`: one-dimensional Euler solver

- constant-`gamma` ideal-gas EOS;
- PCM, componentwise primitive PLM, and time-traced characteristic PLM;
- order-2 or PeleC-style five-point order-4 limited slopes;
- optional pressure/velocity shock flattening;
- Rusanov or qualified single-species PeleC-style Riemann fluxes;
- outflow and periodic boundaries;
- Sod, Shu-Osher, and symmetric planar Sedov-type regressions.

### `pelef2d`: two-dimensional Euler scaffold

- uniform periodic Cartesian mesh;
- x/y directional Riemann fluxes through explicit momentum rotation;
- limited primitive slopes and normal characteristic tracing;
- CTU-style transverse half-step corrections with positivity scaling;
- one unsplit conservative update;
- periodic isentropic-vortex analytical and convergence regressions.

### `pelef3d`: three-dimensional periodic Euler solver

- uniform periodic Cartesian mesh with x/y/z face fluxes;
- Rusanov or qualified single-species PeleC-style Riemann fluxes;
- summed directional spectral-radius CFL selection;
- conservative PCM spatial divergence and transactional SSPRK2 stepping;
- exact x/y/z reduction to the one-dimensional PCM solver;
- periodic entropy-wave convergence, conservation, and deterministic CSV gate.

This baseline 3D application is constant-`gamma` and single-level. Use
`pelef_reactive_3d` for the qualified general-EOS chemistry and transport
path. Molecular transport remains outside this constant-`gamma` baseline;
AMR, restart, MPI decomposition, and EB geometry remain outside both 3D
application boundaries.

### `pelef_reactive_3d`: general-EOS multispecies 3D Euler solver

- elementary seven-species or full ten-species H2/O2 NASA7 thermodynamics;
- periodic x/y/z general-EOS Rusanov, HLLC, or PeleC-style fluxes;
- summed directional CFL selection and transactional SSPRK2 temperature recovery;
- selectable PCM or frozen-composition characteristic PLM with MC/minmod limiting;
- conservation of mass, three momenta, energy, and every species density;
- exact x/y/z reduction and diagonal entropy-wave convergence gates;
- optional elementary/full-H2O2 cell-local chemistry;
- optional viscosity, conduction, mixture-averaged species diffusion,
  barodiffusion, correction velocity, and species-enthalpy transport;
- minimum hydro/parabolic timestep selection and transactional
  `R-T-H-T-R` composition;
- deterministic x-fastest CSV with primitive, `Y_k`, and `rhoY_k` fields.

The executable remains serial, periodic, uniform-grid, and single-level. Its
hotspot cases are internal coupling and conservation regressions, not external
ignition or transport validation.

### Three-dimensional planar embedded-boundary kernel

- exact cell volume fractions, normalized cell and open-face centroids,
  Cartesian face apertures, EB area/centroid, and solid-to-fluid normal
  integrals for one x-, y-, or z-normal plane;
- validator-enforced cell classification, metric bounds, and cellwise
  discrete aperture-divergence/EB-normal identity;
- general-EOS impermeable slip-wall pressure flux and aperture-weighted
  conservative divergence;
- PCM Riemann faces with zero-gradient physical-domain exterior states;
- transactional forward-Euler state and temperature publication, with
  covered cells unchanged;
- conservative FluxRedist through all open x/y/z neighbor faces, with exact
  uniform-RHS preservation and fluid-volume-weighted component conservation;
- a transactional redistributed Euler path that accepts a CFL-0.5 full-grid
  step rejected by the raw volume-fraction-scaled update at `kappa=0.05`;
- zeroth-order weighted StateRedist through the axis-normal regular receiver,
  with target-volume weighting, uniform-state identity, component
  conservation, and transactional EOS recovery;
- one shared face/divergence residual builder feeding raw, FluxRedist, and
  StateRedist Euler paths under identical PCM wall-flux arithmetic;
- an active-cell full-grid CFL selector and installed
  `pelef_reactive_eb_3d` multi-step application that requires stabilization;
- optional masked elementary chemistry composed transactionally as
  `R(dt/2)-H(dt)-R(dt/2)`, with exact covered-cell identity and H/O/N gates;
- optional StateRedist molecular transport with an active-state diffusive
  bound, aperture-weighted internal face fluxes, zero outer/embedded-wall
  transport flux, a species-outflow limiter, transactional SSPRK2 stages, and
  full `R-T-H-T-R` composition;
- versioned single-level checkpoint/restart with complete thermo, reaction,
  transport, geometry, state, and cumulative-diagnostic records plus
  transactional incompatible/truncated-input rollback;
- x-fastest EB CSV with cell classification, volume fraction, physical fluid
  centroid, primitive fields, and every species density;
- exact stationary/tangential uniform-state invariance for x/y/z walls and
  roundoff fluid-volume conservation for a nonuniform case.

The fluid occupies the coordinate-positive side of one non-grid-aligned
plane. An interior plane exactly on a Cartesian face is rejected because this
first data model requires a positive-volume cut cell to own the EB metric.
Curved/oblique geometry, higher-order face reconstruction or StateRedist,
general wall or FluxRedist transport, AMR/reflux, MPI ownership, and a
public 3D EB parallel application remain later increments. Its optional
chemistry is limited to
active cells in the frozen elementary H2/O2/N2 density-sheet problem. The
current public executable is serial, single-level, PCM, forward-Euler hydro,
and exact-axis-plane only. Optional transport is SSPRK2 and StateRedist-only.
It selects the step before the first reaction stage from the minimum hydro and
transport bounds; no general reaction-aware stability guarantee is claimed.

Run both stabilization variants and the independent output audit with:

```bash
./build/release/pelef_reactive_eb_3d cases/reactive_eb_3d/density_sheet_state.nml
./build/release/pelef_reactive_eb_3d cases/reactive_eb_3d/density_sheet_flux.nml
python3 tools/check_reactive_eb_3d.py \
  --state reactive_eb_3d_state.csv --flux reactive_eb_3d_flux.csv

python3 tools/run_and_capture.py \
  --output reactive_eb_3d_chemistry_inert.log -- \
  ./build/release/pelef_reactive_eb_3d cases/reactive_eb_3d/chemistry_inert.nml
python3 tools/run_and_capture.py \
  --output reactive_eb_3d_chemistry_reactive.log -- \
  ./build/release/pelef_reactive_eb_3d cases/reactive_eb_3d/chemistry_reactive.nml
python3 tools/check_reactive_eb_chemistry_3d.py \
  --baseline reactive_eb_3d_state.csv \
  --inert reactive_eb_3d_chemistry_inert.csv \
  --reactive reactive_eb_3d_chemistry_reactive.csv \
  --inert-log reactive_eb_3d_chemistry_inert.log \
  --reactive-log reactive_eb_3d_chemistry_reactive.log

python3 tools/run_and_capture.py --output control.log -- \
  ./build/release/pelef_reactive_eb_3d \
  cases/reactive_eb_3d/transport_control.nml
python3 tools/run_and_capture.py --output transport.log -- \
  ./build/release/pelef_reactive_eb_3d \
  cases/reactive_eb_3d/transport_only.nml
python3 tools/run_and_capture.py --output coupled.log -- \
  ./build/release/pelef_reactive_eb_3d \
  cases/reactive_eb_3d/transport_chemistry.nml
python3 tools/check_reactive_eb_transport_3d.py \
  --control reactive_eb_3d_transport_control.csv \
  --transport reactive_eb_3d_transport_only.csv \
  --coupled reactive_eb_3d_transport_chemistry.csv \
  --control-log control.log --transport-log transport.log \
  --coupled-log coupled.log

./build/release/pelef_reactive_eb_3d \
  cases/reactive_eb_3d_restart/reference.nml
./build/release/pelef_reactive_eb_3d \
  cases/reactive_eb_3d_restart/checkpoint_stop.nml
./build/release/pelef_reactive_eb_3d \
  cases/reactive_eb_3d_restart/restart.nml
ctest --test-dir build/release --output-on-failure \
  -R '^regression_reactive_eb_restart_3d_'
```

### `pelef_amr_reactive_3d`: static two-level reactive 3D AMR solver

- one configurable, strictly interior Cartesian fine patch;
- integer refinement ratio with `r` hydro and `r^2` transport fine subcycles;
- linearly time-interpolated coarse states on all fine-patch boundaries;
- SSPRK2 time-averaged coarse and fine x/y/z face fluxes;
- selectable PCM or characteristic PLM on both levels, with limited linear
  coarse-to-fine ghost interpolation for PLM;
- area/time-averaged six-face reflux and conservative average-down;
- optional viscosity, thermal conduction, mixture-averaged species diffusion,
  barodiffusion (only with species diffusion), correction velocity, and
  species-enthalpy transport;
- transport-specific limited linear coarse-to-fine ghosts and a common
  six-face interface limiter with a conservative fine-side flux correction;
- transactional `R-T-H-T-R` composition when chemistry and transport are
  enabled together;
- composite all-component conservation and transactional whole-step rollback;
- separate deterministic coarse/fine CSV output and an independent checker;
- when transport is disabled, versioned two-level checkpoint/restart with
  configuration and species
  fingerprints, transactional incompatible/truncated-input rejection, and
  exact uninterrupted/restarted output parity; schema 2 fingerprints the
  reconstruction and limiter as well as the preceding solver metadata.
- when transport is enabled, exclusive schema 4 (fixed public `full_h2o2`) or
  schema 5 (selected chemistry) records the ordered gas-transport database,
  controls, cumulative diagnostics, operator identity, and accepted-step
  phase, with transactional read publication.

This serial 3D AMR boundary is periodic, static, and single-patch.
`pelef_mpi_amr_reactive_3d` distributes the same PCM or characteristic-PLM
hierarchy over parent-aligned rank-local x slabs. Root-formatted compatibility
I/O remains rank-neutral and supports transport-disabled schema 2/3 and
transport-enabled schema 4/5 without storing rank ownership. Fixed
elementary/full-H2/O2 chemistry and configure-time selected mechanisms are
available through the same transactional hierarchy split, including molecular
transport and exact changed-rank transport restart. Boundary-touching patches,
regridding, crash-atomic replacement, scalable I/O, CTU/PPM, and EB remain
separate workstream increments.

### `pelef_mpi_amr_reactive_3d`: rank-local sparse 3D AMR verification

- uneven coarse x slabs and parent-aligned fine ownership, including ranks
  with zero fine cells;
- rank-local conserved state and temperature storage throughout CFL, hydro,
  subcycling, reflux, average-down, and temperature recovery;
- two-plane periodic coarse halo propagation and active-neighbor fine halo
  exchange for PCM or characteristic PLM with MC/minmod limiting;
- time-interpolated, limited coarse-to-fine ghost construction and owner-only
  six-face flux-register accumulation;
- fixed elementary/full-H2/O2 `R(dt/2)-H(dt)-R(dt/2)` with source-phase
  average-down, covered-coarse temperature recovery, and zero-fine-rank
  collective participation;
- owner-local SSPRK2 molecular transport with `r^2` fine subcycling,
  time-averaged six-face reflux, common coarse/fine interface limiting, and
  transactional `R-T-H-T-R` full-physics composition;
- exact NASA7, ownership, patch, timestep, solver, reconstruction, and limiter
  consensus before data-dependent collectives, extended to complete ordered
  reaction and source-integrator context when chemistry is enabled;
- root-only schema-2/schema-3 compatibility and schema-4/schema-5 transport
  checkpoint gather/scatter, with no rank count or ownership stored in the
  checkpoint;
- serial/MPI byte parity at 1, 2, 4, and 8 ranks and two-to-four/eight-rank
  restart parity for both PCM and characteristic PLM.

This is a rank-local stepping and memory boundary for one static, periodic
two-level hierarchy. Root-formatted checkpoint/output is deliberately not
described as scalable I/O, and formatted replacement does not claim
crash-atomic persistence. Dynamic topology, physical boundaries, CTU/PPM, 3D EB
AMR, scalable output, and performance scaling remain open.

### `pelef_ms`: passive multispecies Euler solver

- runtime species count with conserved `rho*Y_k` components;
- checked mass-fraction conversion, positivity, and species closure;
- species fluxes tied exactly to the shared mass flux;
- one-dimensional characteristic tracing and two-dimensional CTU transport;
- MultiSpecSod and periodic species-wave regressions.

This older passive path intentionally retains the constant-`gamma` hydro baseline.

### `pelef0d`: thermodynamics and toy-reactor verification

- species molecular weights and NASA7 thermodynamic polynomials;
- mass-based mixture molecular weight, gas constant, `cp`, `cv`, `gamma`, enthalpy, internal energy, and frozen sound speed;
- ideal-gas pressure/density conversion;
- bracketed Newton/bisection inversion from specific internal energy to temperature;
- a synthetic constant-volume two-species isomerization reactor;
- isothermal analytical and adiabatic energy-conservation gates.

### `pelef0d_h2o2`: generated elementary H2/O2 kinetics

- runtime elementary-reaction records with arbitrary reactant/product stoichiometry;
- reversible Arrhenius rates and NASA7 equilibrium constants;
- molar concentrations, progress rates, production rates, and mass-fraction source terms;
- JSON-to-Fortran mechanism generation with a committed-source cleanliness gate;
- an adaptive explicit RK4 constant-volume, adiabatic reactor;
- a seven-species, four-reaction H2/O2/N2 subset selected from Cantera `h2o2.yaml`;
- live trajectory and exact-state production-rate comparison against Cantera 3.2.


### `pelef_transport_probe`: dilute-gas transport verification

- Lennard-Jones transport records for H2, H, O, O2, OH, H2O, and N2;
- Chapman--Enskog pure viscosities and binary diffusion coefficients;
- Wilke mixture viscosity;
- modified-Eucken pure conductivity and Mathur mixture conductivity;
- mixture-averaged species diffusion coefficients;
- live qualification against Cantera 3.2 at four temperature, pressure, and
  composition states.

This is a deliberately qualified ideal-gas subset. It does not yet reproduce
PelePhysics polynomial transport fits, polar corrections, Soret diffusion,
multicomponent diffusion, or bulk viscosity.

### `pelef_reactive_1d`: general-EOS reactive Euler solver

- conserved state `(rho, rho*u, rho*v, rho*w, rho*E, rho*Y_k)`;
- NASA7 composition-dependent pressure, temperature, heat capacities, ratio of specific heats, and frozen sound speed;
- safeguarded conserved-to-primitive recovery through `e(Y,T) -> T` inversion;
- PCM, frozen-composition characteristic PLM, monotone primitive PPM, or
  time-traced frozen-composition characteristic PPM reconstruction;
- selectable Rusanov, general-EOS HLLC, or NASA7 PeleC-style acoustic flux
  with exact species-flux closure;
- SSPRK3 time integration for the semidiscrete primitive-PPM path;
- PeleC-style parabolic profile integration over the `u-c`, `u`, and `u+c`
  waves for characteristic PPM;
- optional PeleC one-dimensional shock flattening and bounded
  Colella--Woodward contact steepening on characteristic PPM;
- periodic or outflow boundaries;
- cell-local adiabatic constant-volume chemistry;
- optional shear viscosity, Fourier heat conduction, mixture-averaged species
  diffusion, barodiffusion, correction velocity, and species-enthalpy flux;
- explicit SSPRK2 transport with a parabolic timestep gate;
- symmetric reaction--transport--hydro--transport--reaction splitting;
- homogeneous-reactor reduction, smooth density/composition-wave convergence,
  discontinuous material-contact resolution, primitive/characteristic PPM
  convergence, periodic strong-shock flattening, and nonuniform
  reactive-hotspot regressions.


### `pelef_reactive_2d`: general-EOS reactive CTU solver

- uniform periodic Cartesian mesh with the same conserved reactive state as the 1D path;
- composition-dependent NASA7 pressure, temperature, heat capacities, and frozen sound speed;
- directional general-EOS Rusanov, HLLC, or PeleC-style acoustic fluxes
  through explicit momentum rotation;
- PCM, frozen-composition characteristic PLM, or time-traced
  frozen-composition characteristic PPM in both coordinate directions;
- optional bounded contact steepening and PeleC-style shock flattening on the
  characteristic-PPM normal predictor;
- provisional face fluxes, conservative CTU transverse half-step corrections, and EOS-based positivity scaling;
- species, momentum, and total energy corrected together so `sum(rho*Y_k)=rho` remains coupled to the hydro update;
- cell-local chemistry with symmetric reaction--transport--hydro--transport--reaction splitting;
- optional x/y Newtonian viscosity, Fourier conduction, mixture-averaged species diffusion, barodiffusion, correction velocity, and species-enthalpy flux;
- exact diagonal density/composition-wave convergence, x/y one-dimensional
  reduction, material-contact sharpening, oblique strong-shock flattening,
  periodic vortex, and reacting-hotspot regressions.

The reactive 2D path supports matched periodic pairs, slip/no-slip walls,
adiabatic/isothermal wall temperatures, fixed-state inflow, and zero-gradient
outflow. Walls are species-impermeable by default or may impose a prescribed
zero-net-mass species-conversion flux with consistent species-enthalpy
transport. Periodic cases retain the qualified characteristic-PPM plus CTU
path; physical faces use boundary-aware ghost reconstruction, an exact
impermeable inviscid wall flux, and boundary-aware molecular transport.

The reactive applications can select either the verified seven-species,
four-reaction elementary subset or the full ten-species, 29-reaction H2/O2
mechanism with third-body, falloff, Troe, and adaptive implicit chemistry. The
characteristic projection is a qualified frozen-composition ideal-gas-mixture
approximation, not full PeleC/PelePhysics general-EOS characteristic parity.

### Embedded-boundary geometry foundation

`eb_geometry_2d_mod` converts a nodal level set into bounded Cartesian
cell-volume fractions, x/y face-area fractions, and regular/cut/covered cell
types. Positive level-set values define fluid. Each cell uses two affine
triangles, so planar interfaces are exact and curved interfaces converge under
refinement. Cut cells also carry the physical embedded-boundary length,
centroid, and a unit normal directed from solid toward fluid. The reactive EB
kernel recovers the general-EOS pressure and applies the
stationary impermeable slip-wall momentum flux on an arbitrarily oriented cut
wall. It also converts the integrated wall flux into a volume-normalized source
for each cut cell without mass, energy, or species leakage. Curved-cell force
uses the integrated interface-normal vector rather than multiplying total
length by an averaged unit normal. The conservative divergence combines this
wall contribution with shared Cartesian fluxes weighted by open-face fraction,
and preserves a uniform stationary pressure field for planar and circular
embedded boundaries. Molecular transport uses the physical distance from the
fluid-volume centroid to the embedded-wall centroid along the solid-to-fluid
normal. Isothermal walls contribute Fourier heat flux; no-slip walls contribute
the corresponding normal Newtonian traction and viscous work at the configured
wall velocity. Slip, adiabatic, and every species-impermeable default remains
exactly zero. Small-cell time integration now has a conservative
first-order FluxRedist path: it blends a cut-cell update with its
volume-weighted face-connected neighborhood, redistributes the removed
extensive update, and commits a forward update only after every active reactive
state passes EOS recovery. Weighted StateRedist forms normal-directed
neighborhoods with target volume fraction `0.5`, accounts for cells shared by
overlapping neighborhoods, conserves every volume-weighted state component,
and applies the same transactional EOS gate. Its default `max_order=0` path
retains neighborhood averages. Selectable `max_order=2` stores normalized
fluid-volume centroids, fits limited linear `Qhat` slopes over active
neighborhoods, evaluates them at every merge recipient, and adds a
conservation-compatible recipient maximum-principle limiter.
A selectable PCM or frozen-composition characteristic-PLM EB hydro path now
constructs reactive Riemann fluxes only on open Cartesian faces, uses
zero-gradient domain faces, linearly interpolates face-center fluxes to the
open-face centroid, combines them with the integrated slip-wall pressure
force, and completes the step through weighted StateRedist and EOS recovery.
PLM slopes use only two-sided active-cell stencils; cells beside covered or
outer cells fall back locally to zero slope. The `pelef_reactive_eb_2d`
application now reads plane or circular geometry from a namelist, initializes a
general-EOS multispecies state, advances to a requested final time with an
active-cell CFL limit, reports volume-weighted diagnostics, and writes cell
geometry and primitive fields to CSV. This qualified runnable path has
zero-gradient outer faces and optional active-cell chemistry. Its
Strang sequence applies half reactions only to active cells and can insert
SSPRK2 molecular-transport half steps around the EB hydro transaction while
leaving covered cells bitwise unchanged. Transverse reconstruction settings
remain rejected instead of silently ignored.

The EB and AMR foundations now meet at one qualified static two-level transfer.
An aligned rectangular fine patch is restricted with fine fluid-volume weights,
the corresponding composite integral counts uncovered coarse cells and fine
cells exactly once, and reactive restriction recovers every active parent
temperature before committing either state or temperature. Covered parents
retain their original reactive data. This is a serial synchronization kernel,
not yet a time-advancing EB AMR application.

The same static hierarchy now owns an EB-aware flux register. Coarse and fine
steps accumulate time-integrated face fluxes independently; open-face fractions
and physical subface lengths produce one correction on each exterior coarse
cell. A regular cell receives that correction directly. A cut cell keeps its
fluid-volume share and redistributes the remainder over connected 3-by-3
neighbors, with any share landing below the fine patch transferred to all of
that parent's fine children. Reactive re-reflux commits both levels and both
temperature fields only after every active cell passes EOS recovery.

The static hierarchy can now initialize its fine rectangle by piecewise-
constant coarse injection and advance reactive EB hydrodynamics over one
coarse interval. The coarse level takes one step while the fine level takes
`r` steps. Each fine-patch boundary face obtains a time-interpolated conserved
state from the adjacent coarse cell, followed by EOS temperature recovery.
The fluxes that actually advanced both levels feed the EB flux register;
re-reflux and reactive average-down then synchronize the hierarchy in one
transaction.

`pelef_reactive_eb_amr_2d` makes that hierarchy runnable from one input file.
It builds a strictly internal fine rectangle, initializes it from the coarse
state, selects each coarse timestep from both level CFL limits, advances until
the clipped final time, and writes separate synchronized coarse and fine
geometry/state CSV files. Optional solution-driven regridding tags active
coarse cells by relative and absolute temperature jumps, buffers their bounding
rectangle, averages the old fine patch down, injects the new patch from coarse
data, and retains every overlapping same-resolution fine cell exactly.
When configured to remove an untagged patch, the driver conservatively collapses
the child into the root, releases its arrays and geometry, advances with the
single-level CFL and hydro path, and re-creates the patch by PCM when tags
return. Fine CSV output is omitted while the child is inactive.

With chemistry enabled, each active AMR interval applies a reaction half-step
on both levels, the existing subcycled EB hydro/reflux transaction, a second
reaction half-step, and final fine-to-coarse average-down. Covered cells are
masked from chemistry. A root-only lifecycle interval reuses the qualified
single-level EB Strang path, and any chemistry, hydro, or EOS failure leaves the
complete hierarchy unchanged.

The EB AMR input can request periodic checkpoints, stop immediately after a
successful write, or restart from a prior file. A checkpoint records the
coarse state, optional fine state and current patch bounds, time/step/regrid
metadata, minimum accepted timestep, base density, species ordering, and the
physics and regrid settings required for compatible continuation. Restart
rebuilds both EB geometries from the input, recovers active-cell temperatures
from the conserved state, and publishes the hierarchy only after the complete
file and end marker pass validation. Final time, maximum steps, output paths,
and checkpoint controls may change between runs.

The two-level EB kernel can also cluster disconnected tag components into a
deterministically ordered collection of separated fine rectangles. Buffer and
minimum-size expansion are applied per component, nearby candidates are
coalesced to preserve the 3-by-3 redistribution neighborhood contract, and
clamped to the complete root domain. Boundary cells use available one-sided
temperature differences, so a single- or multipatch plan may meet any outflow
side. A patch-set transaction
supports PCM creation, exact old/new fine-overlap retention, conservative
average-down, composite integration, movement, repartition, and removal.
Hydrodynamics advances the root once, subcycles every child with its own EB
flux register, then refluxes and averages down all children. The matching
Strang transaction applies active-cell chemistry on the root and every child
around that hydrodynamic interval and restores the complete hierarchy after
any rejected stage.

Setting `multipatch_enabled = .true.` connects that patch set to the public
`pelef_reactive_eb_amr_2d` lifecycle. The input exposes the maximum tag gap,
the initial and periodic planners can replace the configured seed rectangle
with zero or more children, the timestep is the minimum root/child CFL limit,
and accepted steps retain the existing regrid cadence and counters. Fine CSV
paths receive deterministic `_patch0001`, `_patch0002`, ... suffixes.

Multipatch mode has a separate versioned formatted checkpoint schema. It stores
the root, ordered child count, every child's actual bounds and state, accepted
time and timestep, step/regrid counters, base density, species order, and a
strict physics/topology compatibility signature. Restart rebuilds every EB
geometry, recovers active temperatures through the EOS, validates the complete
set and end marker in private candidates, and only then publishes the restored
hierarchy. The earlier single-patch schema and its inputs remain unchanged.

Setting `three_level_enabled = .true.` instead constructs one static middle
rectangle from root indices and one finest rectangle from middle indices. The
finest rectangle retains a two-cell middle margin. The public timestep is the
minimum root-equivalent stability limit from all three levels, and every
accepted interval uses recursive subcycling, independent interface registers,
the EB-cut conservation closure, active-cell Strang chemistry, and final
deepest-first synchronization. Successful completion writes distinct root,
middle, and finest CSV files. Scheduled and final three-level checkpoints use
a dedicated schema and may stop and resume the same hierarchy without changing
the established single-patch or patch-set formats. Three-level mode remains
mutually exclusive with multipatch siblings. Its dynamic path always retains
the finest patch and ignores finest-level tags outside the two-cell-safe
planning region. By default the middle patch remains fixed. Setting
`dynamic_parent_regridding = .true.` also replans that patch from root tags,
rebuilds the complete nested hierarchy atomically, and uses dynamic checkpoint
schema 4 to restart both actual patch descriptors and the regrid cadence.

Unsplit transverse prediction, fourth-order StateRedist slopes,
periodic/ghost-cell neighborhoods, catalytic embedded-wall species transfer,
higher-order wall-normal gradients, public selection of limited-linear
coarse-to-fine initialization, and multiple dynamic parents at one fixed-depth
level are not yet connected. The fixed-depth public EB
AMR application remains serial and owns either restartable sibling rectangles
or an explicit three-level hierarchy with optionally dynamic parent and finest
patches;
the separate arbitrary-depth application provides the qualified sparse-MPI
lifecycle.

The separate EB patch-tree core now owns reactive conserved state and
temperature on arbitrary-depth, branching topology. It initializes children
from their actual parents, synchronizes deepest-first, evaluates the complete
composite conserved vector, and transactionally migrates fields through a
whole-tree topology replacement. Same-resolution physical overlap is retained
only after local EB geometry checks; EOS or conservation failure leaves the
accepted tree unchanged. It also evaluates the active-cell CFL limit on every
runtime node containing fluid and scales each local limit by the cumulative
refinement product before publishing one stable root interval. Reactive EB
hydrodynamics now walks the same runtime tree recursively: every child takes its
exact ratio subcycles, owns one flux register, refluxes into its actual parent,
and participates in subtree conservation closure and deepest-first final
synchronization. Active-cell chemistry advances every node and composes two
reaction half-steps around that recursive hydro transaction with one atomic
tree commit. SSPRK2 molecular transport also recurses over the runtime tree,
subcycles each relation, refluxes diffusive fluxes, closes every refined
subtree, and blends both Euler stages before commit. One combined `R-T-H-T-R`
entrypoint now applies both reaction halves, both SSPRK2 transport halves, and
the recursive hydro interval on the same private tree candidate. A public
target-time loop recomputes the all-node hydro/transport stability limit before
every step, clips exactly to the requested stop time, and commits the tree,
clock, step count, limiter minimum, and per-level physics counts together.
The same serial tree now synchronizes accepted fields, tags temperature
gradients independently on every prospective parent, clusters disconnected
features, constructs caller-defined EB child geometry, and transactionally
rebuilds or collapses the arbitrary-depth topology. The same tree writes a
single composite CSV in which cells covered by finer AMR patches are omitted.

An MPI ownership descriptor maps every arbitrary-depth tree node to a
deterministic rank using allocated cells and optional subcycle-weighted work.
The sparse numerical representation allocates state and temperature only on
that node's owner while retaining replicated topology and ownership metadata.
Explicit materialization reconstructs a complete tree when required, and a
new owner map migrates changed nodes directly from old owner to new owner into
a private candidate. Combined hydro and explicit-transport stable-step
selection now evaluates only owner-local active nodes, converts every result
to the root interval, and reduces the global minimum without materialization.
Chemistry also advances only owner-local nodes, then synchronizes the hierarchy
deepest-first by sending child state directly to a distinct parent owner before
average-down. Whole-tree and selected-subtree composite conserved integrals
now recurse over the replicated topology, exclude refined parent cells, sum
only owner-local node contributions, and reduce one conserved vector without
materializing nonowned fields. Recursive hydro now advances that same sparse
tree on owners, routes compact parent-time context and fine fluxes across
distinct-owner edges, applies ordered reflux/average-down, and closes each
refined subtree before one atomic commit. Recursive SSPRK2 molecular transport
now uses the same direct owner route for both Euler stages, performs its blend
and EOS temperature recovery only on owners, and synchronizes deepest-first
without materialization. One outer sparse transaction now composes owner-local
chemistry half-steps, SSPRK2 transport half-steps, and recursive hydro as
`R-T-H-T-R`, committing fields, limiter, advances, and transfer counts only
after every stage succeeds. A public sparse clock now recomputes the owner-
local stable interval before every step, clips the last step exactly to the
requested target time, and commits fields, time, step count, minima, advances,
and transfers together. Owner-local MPI arbitrary-depth temperature tagging
now reduces compact per-parent rectangles, rebuilds EB geometry on every rank,
recomputes deterministic ownership, and migrates retained overlap directly
between old and new owners without materializing a complete numerical tree.
The serial arbitrary-depth tree now writes and reads a self-describing
formatted checkpoint containing its complete topology, EB metrics, fields,
and lifecycle metadata. Sparse MPI writes now gather numerical nodes only to a
selected root, while restart broadcasts geometry, recomputes ownership for the
current rank count, and scatters fields directly to their new owners. Sparse
MPI composite output reuses the selected-root gather boundary, invokes the
serial writer only there, and leaves complete numerical output unallocated on
all other ranks.

The public `pelef_reactive_eb_patch_tree_2d` executable connects this serial
tree to the existing `&reactive_2d`, `&embedded_boundary`, and `&eb_amr`
inputs. `patch_tree_maximum_levels` bounds recursive temperature-tag planning;
the established refinement ratio, clustering controls, physics switches,
checkpoint paths, and output path drive the same qualified core operations.
The application starts from a root-only EB field or the self-describing tree
checkpoint, optionally regrids at initialization and at the configured step
cadence, advances one stable root interval transactionally, and publishes one
composite CSV at completion or checkpoint stop.

The public patch-tree lifecycle is also qualified across a real application
checkpoint boundary. A four-level dynamic reference run is compared with a
run stopped after its first scheduled checkpoint and a continuation loaded by
a separate process. The restart restores the arbitrary-depth geometry, state,
time, root-step and regrid counters, minimum accepted timestep, minimum
transport limiter, and original composite-integral baseline; its final
composite topology, numerical fields, and cumulative diagnostics match the
uninterrupted result.

The installed `pelef_mpi_reactive_eb_patch_tree_2d` executable exposes the
same input-driven lifecycle through sparse MPI ownership. Recursive tagging,
timestep selection, full physics, regridding, checkpoint/restart, integrals,
and composite output operate on owner-local node fields. The workload exponent
is configurable from `&eb_amr`, and 1, 2, 4, and 8 ranks produce the same
four-level composite result.

The replicated MPI-owner EB AMR hydro path now decomposes the root update over
its distributed y-tiles. Each tile owner advances a bounded six-row halo band,
publishes only its owned cells and uniquely assigned faces, and participates in
a collective root assembly. The previous single root-owner full-level advance
and four full-root broadcasts are absent from this path. Fine children retain
their established owner subcycling and deterministic reflux order. The sparse
root path remains a later conversion boundary.

### MPI one-dimensional verification

With `PELEF_ENABLE_MPI=ON`, five domain-decomposed executables verify:

- uneven non-replicated block decomposition and periodic halo exchange;
- conservative multispecies Euler transport;
- distributed adaptive implicit full-H2/O2 chemistry scheduling;
- general-EOS molecular transport with viscosity, conduction, barodiffusion,
  correction velocity, and species-enthalpy transport;
- transactional reaction--transport--hydro--transport--reaction splitting;
- ordered gather output, global timestep/conservation reductions, and
  complete-field parity for 1, 2, 4, and 8 ranks.

A sixth MPI executable exercises sparse AMR distribution. Compact hierarchy
and owner metadata are replicated, but each root/fine field payload exists
only on its deterministic work-weighted owner. Patch work may use raw cell
count, hyperbolic `r` subcycling, or parabolic `r^2` subcycling; explicit and
tag-driven regrids preserve the selected model. Rank-local hyperbolic and
parabolic stability limits are reduced to one communicator-wide coarse-step
limit without gathering patch fields. Point-to-point transfers cover
same-level halos, parent/child ghost data, boundary fluxes, shared-face
corrections, average-down, explicit regrid prolongation, retained overlap, and
owner migration. Chemistry, recursive hydro, parabolic `r^2` transport, and
the transactional `R-T-H-T-R` interval run on owners alone. Both explicit and
solution-tagged topology rebuilds remain field-sparse; owner-local tagging
shares only compact integer plan metadata. The 1/2/4/8-rank gates compare the
gathered hierarchy, fields, ghosts, counters, conservation, and rollback with
the serial patch-tree implementation.

The seventh executable, `pelef_mpi_amr_reactive_1d`, is the public sparse AMR
driver. It reads the reactive namelist, builds owner-local tagged patch trees,
selects distributed hydro/transport timesteps, advances `R-T-H-T-R`, regrids
at the requested cadence, and writes an ordered composite AMR CSV. Persistent
field payloads remain globally single-copy during the time loop; a complete
tree is materialized for final diagnostics/output and scheduled checkpoints.
Checkpoints store no owner map: restart rebuilds deterministic ownership for
the active communicator, allowing a two-rank run to resume on four or eight
ranks.

The eighth executable, `pelef_mpi_eb_amr_patch_2d`, establishes the first
two-dimensional EB AMR distribution boundary. It partitions the root into
contiguous y-tiles, assigns separated fine siblings as independent entities,
and balances raw, hyperbolic, or parabolic subcycle-weighted work with 64-bit
accounting. Geometry and patch metadata remain replicated, while root-tile and
child state/temperature payloads are synchronized from one authoritative
owner. Its 1/2/4/8-rank gate checks collective topology agreement, exact
ownership accounting, authoritative payload recovery, and transactional
rejection of invalid or inconsistent work models. The same gate now advances
active-cell chemistry exactly once on each root-tile or child owner, reaches
serial patch-set parity after owner broadcasts and fine-to-root average-down,
and rolls every rank back after a late owner-side reactor rejection. Reactive
EB hydrodynamics now advances the complete root level once on its exclusive
physics owner and advances every child on its patch owner with ratio
subcycling, owner-local flux registers, sequential reflux, and final
fine-to-root average-down. Its result matches the serial multipatch path and a
late child-owner failure leaves every rank unchanged.

### One-dimensional AMR

The AMR layer provides:

- a static two-level hierarchy with an interior refined patch;
- an arbitrary-depth hierarchy foundation built from nested level interfaces;
- arbitrary component counts and integer refinement ratios;
- MC-limited, conservative piecewise-linear prolongation;
- volume-average restriction and covered-cell synchronization;
- refinement-ratio level subcycling;
- time-integrated coarse/fine flux registers and reflux;
- a composite-integral gate proving conservation across both patch interfaces;
- component-selectable normalized-gradient tagging with absolute noise floors;
- buffered, minimum-width single-patch planning with explicit boundary checks;
- conservative patch creation, movement, resizing, and removal;
- exact retention of same-resolution fine data where old and new patches overlap;
- reactive general-EOS states advanced on both coarse and fine levels;
- refinement-ratio fine-level subcycling with time-interpolated coarse ghosts;
- coarse/fine flux accumulation, reflux, and average-down every coarse step;
- symmetric chemistry--hydro--chemistry splitting across the AMR hierarchy;
- transactional hierarchy rollback and periodic solution-driven regridding;
- ordered composite CSV output with exact domain-coverage checks;
- optional limited primitive-variable PLM with SSPRK2 time advancement;
- time-averaged SSPRK2 interface fluxes used consistently for reflux;
- four-layer physical and coarse/fine ghost storage for PPM stencils;
- primitive or time-traced characteristic PPM with SSPRK3 level advancement;
- optional PeleC-style hybrid WENO5-JS, WENO5-Z, WENO7-Z, or WENO3-Z edge
  reconstruction inside the characteristic-PPM path;
- SSPRK3 effective interface fluxes used consistently for reflux;
- conservative, MC-limited parent-to-fine ghost interpolation at subcycle
  midpoint times;
- AMR viscosity, Fourier conduction, and mixture-averaged species diffusion;
- parabolic fine subcycling with time-interpolated coarse transport ghosts;
- diffusive flux-register reflux and covered-cell average-down;
- symmetric reaction--transport--hydro--transport--reaction composition;
- cumulative multilevel subcycle schedules and deepest-to-root synchronization;
- a four-level, mixed-ratio conservation gate across every interface;
- arbitrary-depth reactive state and temperature ownership;
- recursive hydro subcycling and diffusive `r^2` subcycling at every depth;
- recursive chemistry and reaction--transport--hydro splitting;
- a three-level reactive conservation, positivity, closure, and synchronization
  gate;
- a runtime `amr_max_levels` limit with tag-driven nested hierarchy creation;
- periodic conservative multilevel hierarchy rebuilds during simulation;
- recursively ordered, exact-coverage multilevel composite CSV output;
- a runnable three-level hotspot case and structural output gate;
- physical-coordinate overlap transfer across changed multilevel hierarchies;
- exact retention of aligned old fine state and temperature data;
- nested fine patches touching an outflow physical boundary;
- physical-side PPM/WENO ghost fill with reflux restricted to the remaining
  coarse/fine interface;
- configured two-level EB fine patches touching an outflow physical boundary,
  with zero-gradient physical-side exterior state and no physical-side reflux;
- ordered sets of disjoint fine patches over one parent level;
- disconnected-tag clustering with deterministic buffer/minimum-width
  expansion and automatic coalescing of adjacent candidates;
- set-wide conservative prolongation, average-down, reflux, and composite
  integration without double counting covered parent cells;
- conservative patch-set movement, repartition, removal, and exact retention
  of same-resolution fine overlap across old/new patch intersections;
- fixed two-level reactive WENO7-Z hydro subcycling on two separated fine
  patches with per-patch flux registers;
- fixed two-level multipatch reaction--transport--hydro splitting with `r^2`
  diffusive subcycling, per-patch reflux, and transactional rollback;
- tag-driven two-level multipatch creation, movement, repartition, removal,
  overlap retention, runtime statistics, and ordered composite CSV output;
- arbitrary-depth separated patch trees with explicit parent ownership,
  mixed per-level refinement ratios, recursive conservative prolongation,
  deepest-to-root average-down, and exact composite integration;
- parent-owned patch-tree flux registers with transactional deepest-to-root
  reflux and covered-cell synchronization across every branch;
- static arbitrary-depth reactive patch trees with recursive PCM hydro
  subcycling, time-interpolated parent ghosts, per-child flux accumulation,
  reflux, and average-down at every branch;
- a four-level `1/2/3/2`-patch reactive conservation, synchronization,
  positivity, closure, and exact subcycle-count gate;
- symmetric chemistry--hydro--chemistry splitting over every patch-tree node,
  with deepest-to-root post-reaction synchronization and whole-tree rollback;
- a chemistry-on versus hydro-only branch comparison proving reactive species
  evolution while retaining composite mass, momentum, and total energy;
- recursive molecular transport with cumulative `r^2` child subcycling,
  time-interpolated parent ghosts, per-child diffusive flux registers, reflux,
  and average-down at every branch;
- full reaction--transport--hydro--transport--reaction composition with a
  four-level transport-on versus transport-off conservation gate;
- transactional runtime patch-tree rebuilds from supplied branching plans,
  with physical-coordinate overlap transfer across changed parent ownership;
- no-op, moved-tree conservation, deepest exact state/temperature retention,
  counter preservation, and invalid-plan rollback gates;
- per-parent normalized-gradient tagging and deterministic disconnected-tag
  clustering for automatic arbitrary-depth branching plans;
- tag-driven transactional tree rebuilds, including root-only creation,
  maximum-depth branching, unchanged-plan no-op, and invalid-request gates;
- independently owned adjacent children with parent-local same-level exchange
  for face and four-layer PPM/WENO ghost data;
- single owned time-integrated hydro and diffusive fluxes at fine/fine faces,
  with those internal faces excluded from coarse/fine reflux;
- adjacent PPM hydro and molecular-transport conservation, synchronization,
  exact exchange, and subcycle-accounting gates;
- deterministic cell-, hyperbolic-subcycle-, or parabolic-subcycle-weighted
  MPI ownership for every tree patch, collective hierarchy/work-model
  consensus rejection, owner-authoritative patch synchronization, and
  four-layer cross-rank adjacent-sibling halo gates;
- owner-only patch-tree chemistry with one global advance per patch, serial
  field parity, deepest-to-root synchronization, conservation, and global
  transactional rollback gates;
- owner-only recursive MPI patch-tree hydro with exact per-owner subcycle
  accounting, cross-owner adjacent PPM flux reconciliation, serial field
  parity, conservation, and global transactional rollback gates;
- owner-only recursive MPI patch-tree molecular transport with exact
  parabolic subcycle accounting, cross-owner shared diffusive fluxes, serial
  parity, conservation, and global transactional rollback gates;
- a transactional owner-only MPI `R-T-H-T-R` interval with complete
  bookkeeping synchronization, serial full-field parity, exact operator-call
  accounting, missing-database rejection, and rollback after a later-stage
  failure;
- rank-local sparse MPI AMR patch payloads, exact owner scatter/gather,
  one-copy global storage accounting, and same-hierarchy owner migration;
- direct owner-only chemistry on sparse AMR payloads with distributed
  average-down, parent/child ghost fill, adjacent PPM ghost replacement, and
  exact collective rollback;
- targeted point-to-point sparse MPI EB restriction from each child owner only
  to intersecting root tile owners, with exact transfer accounting and
  transactional rollback;
- owner-tiled direct sparse MPI EB hydro with point-to-point six-row halo
  exchange, bounded tile-local work, direct tile-to-child state/flux and
  correction routing, zero root-result/scatter traffic, exact accounting, and
  serial parity;
- owner-tiled direct sparse MPI EB SSPRK2 transport with point-to-point
  six-row halos, seam-isolated finite periodic-edge bands, targeted
  result/scatter and child traffic, zero-traffic tile-local final blending,
  exact accounting, and serial parity;
- patch-local EB flux-register storage and compact sparse MPI transport reflux
  correction round trips with unchanged transactional ordering;
- owner-local sparse MPI EB hydro/transport timestep selection with no root
  field traffic, fine-to-coarse subcycle scaling, communicator-minimum
  reduction, serial timestep parity, and collective rejection;
- public sparse MPI EB multi-step `R-T-H-T-R` advancement with a freshly
  selected stable interval per step, exact final-time clipping, committed-only
  clock and diagnostic publication, and serial full-field parity;
- transactional explicit sparse MPI EB topology rebuilding with serial
  overlap/prolongation parity, direct restriction/PCM/overlap owner traffic,
  deterministic owner recomputation, one-copy post-regrid storage, and
  complete invalid-control rollback;
- scheduled temperature-tagged sparse MPI EB topology rebuilding with
  root-owner planning, compact metadata broadcast, caller-defined geometry,
  direct owner migration, exact transfer accounting, serial dynamic-loop
  parity, and whole-step rollback;
- targeted root-only sparse MPI EB materialization with exact field parity,
  one packed send per remote root tile or child, unallocated non-root outputs,
  and collective invalid-payload rollback;
- writer-root-only sparse MPI EB formatted checkpoint and root/child CSV
  publication, serial checkpoint round-trip compatibility, collective I/O
  status, and exact successful-transfer accounting;
- root-only formatted checkpoint read with direct root-to-owner sparse restart
  scatter from a geometry-only replicated topology descriptor, exact
  per-remote-entity traffic, rank-local field parity, and collective
  metadata/I/O failure rollback;
- direct recursive hydro on sparse AMR payloads with mixed-ratio subcycling,
  replicated flux-register metadata, owner-local reflux/average-down,
  cross-owner PPM face reconciliation, and exact rollback;
- direct recursive molecular transport on sparse AMR payloads with cumulative
  `r²` subcycling, diffusive flux registers, cross-owner shared-face
  reconciliation, and exact rollback;
- direct SSPRK2 molecular transport on sparse MPI EB AMR payloads with
  owner-local child subcycling and reflux, distributed cut-interface
  conservation closure, serial parity, limiter parity, and exact rollback;
- a direct sparse `R-T-H-T-R` transaction with owner-only stage execution,
  exact call accounting, serial parity, missing-database rejection, and outer
  rollback after a later-stage failure;
- an end-to-end sparse MPI EB AMR `R-T-H-T-R` transaction that composes direct
  sparse chemistry, SSPRK2 transport, and hydro without a replicated child
  compatibility window;
- transactional topology-changing sparse regrid with rebuilt owner maps,
  exact overlap retention, one-copy persistent storage, serial parity, and
  invalid-plan rollback;
- tag-driven sparse regrid through four levels with disconnected-feature
  clustering, unchanged-plan no-op behavior, conservation, and invalid-tag
  rollback;
- packed point-to-point same-hierarchy owner migration with one direct message
  per changed patch and exact state, temperature, and ghost reconstruction;
- packed point-to-point adjacent sparse halo exchange with one bidirectional
  payload per cross-owner sibling face and no traffic on unrelated ranks;
- direct child-owner to parent-owner sparse interior transfer for chemistry
  average-down and hydro/transport synchronization;
- direct parent-state fanout only to distinct remote owners that require the
  state for sparse child ghost refresh;
- broadcast-free sparse recursive hydro and transport with owner-local flux
  registers, direct interval/flux/correction payloads, and reduced counters;
- an input-driven sparse MPI AMR driver with configurable subcycle-weighted
  ownership, stop-time clipping, periodic owner-local regridding, final
  conservation diagnostics, and ordered composite patch-tree CSV output;
- a moving-contact gate demonstrating lower AMR error than PCM.

For PCM/PLM, the reactive AMR application retains its overlap-preserving
two-level path when `amr_max_levels = 2`. PPM selects the multilevel engine at
any configured depth so every level owns the wider stencil state. The same
engine provides tag-driven arbitrary-depth state ownership, recursive
advancement, periodic hierarchy rebuilds, and composite output. A changed
multilevel hierarchy is conservatively averaged to the root before nested
patches are rebuilt, then old fine state and temperature data are copied
wherever old/new level spacing and physical cells align. Changed refinement
ratios fall back to conservative prolongation. With
`amr_multipatch_enabled = .true.`, the public application uses a tag-driven
two-level patch set, clusters disconnected tags, periodically rebuilds the
set, retains aligned fine overlap, and writes every uncovered parent or fine
cell exactly once. A separate static patch-tree engine permits each parent
patch to own zero or more separated children at arbitrary depth and advances
those patches recursively with level-ratio hydro subcycling. Every local child
accumulates its own coarse/fine boundary fluxes; reflux and average-down occur
after its fine subcycles, and failures restore the complete tree. The
patch-tree engine can synchronize to the root, tag and cluster every
prospective parent independently, and rebuild the resulting arbitrary-depth
branching plan transactionally. Independently owned adjacent siblings exchange
same-level ghosts, reconcile each shared time-integrated interface flux, and
exclude that internal side from reflux. The first MPI distribution bridge
assigns a unique owner to every patch and communicates authoritative fields
and adjacent halos while retaining replicas on all ranks. Chemistry and the
recursive hydro patch kernel execute only on those owners with collective
acceptance and rollback. Owner-authoritative start-state and face-flux
broadcasts let every replica apply the existing subcycling, shared-flux,
reflux, and average-down rules deterministically. Owner-only molecular
transport uses the same structure with `r^2` subcycling. One outer MPI
transaction now composes chemistry, transport, hydro, transport, and chemistry
and restores its synchronized starting tree after any rejected stage. The
sparse MPI bridge now removes non-owned patch payloads, reconstructs an exact
replica on demand, and moves payload ownership when a same-hierarchy owner map
changes. Chemistry is the first operator to run directly on sparse payloads:
child interiors are synchronized
deepest-to-root for average-down and parent/sibling data are streamed for
ghost refresh without materializing a complete tree. Recursive hydro now uses
the same sparse storage boundary while streaming interval states and fluxes;
molecular transport now follows the same structure with its parabolic `r²`
schedule. Their combined `R-T-H-T-R` transaction now remains sparse from its
outer backup through final acceptance, so normal full-physics advancement no
longer materializes a complete replica. Explicit-plan topology changes now
construct new patches on their owners, stream parent-prolongated interiors,
and transfer fine overlap directly from old owners to new owners. The same
sparse API now derives arbitrary-depth plans on parent owners, replicating only
compact integer topology metadata while building candidate fields on their
owners. Same-hierarchy owner changes now send one
packed patch directly from the old owner to the new owner. Adjacent sparse
siblings likewise exchange only the one- or four-layer boundary payload needed
by their two owners. Child interiors used by average-down and synchronization
now move directly from each child owner to its parent owner. Parent interval
start/end states now reach only distinct child owners, child boundary fluxes
return directly to the parent owner, and shared-face corrections return only to
the affected child owner. Level counters synchronize from owner-local deltas
once per physics stage. Final parent-to-child ghost refresh likewise sends one
parent state to each distinct remote child owner. Sparse physics and both
explicit-plan and tag-driven regrid contain no all-rank field replica. The
public sparse MPI AMR driver composes those APIs into a complete run; only the
initial root state and final diagnostic/output tree are materialized. A
periodic child may
touch a physical boundary only when it covers the full parent domain;
one-sided periodic refinement remains excluded because it crosses the periodic
seam.

## Build and test

Requirements: CMake 3.23 or newer, a Fortran 2018 compiler, and Python 3.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
ctest --test-dir build --output-on-failure
```

MPI verification additionally requires an MPI implementation with the Fortran
2008 `mpi_f08` module. The MPI wrapper must use the same Fortran compiler
family as `CMAKE_Fortran_COMPILER`; configuration fails early when the module
probe is incompatible:

```bash
cmake -S . -B build-mpi \
  -DCMAKE_BUILD_TYPE=Release \
  -DPELEF_ENABLE_MPI=ON
cmake --build build-mpi --parallel
mpiexec -n 4 ./build-mpi/pelef_mpi_reactive_1d mpi_reactive_np4.csv
mpiexec -n 4 ./build-mpi/pelef_mpi_amr_reactive_1d \
  cases/mpi_sparse_amr_hotspot/hotspot.nml sparse_amr_np4.csv
mpiexec -n 4 ./build-mpi/pelef_mpi_eb_amr_patch_2d
mpiexec -n 4 ./build-mpi/pelef_mpi_amr_reactive_3d \
  cases/amr_reactive_3d/amr_entropy_wave.nml mpi_amr_3d_np4
```

To enable the live Cantera reference gate:

```bash
python3 -m pip install cantera==3.2.0
cmake -S . -B build-cantera \
  -DCMAKE_BUILD_TYPE=Release \
  -DPELEF_ENABLE_CANTERA_REFERENCE=ON
cmake --build build-cantera --parallel
ctest --test-dir build-cantera --output-on-failure
```

Regenerate the pinned full-H2/O2 bundle explicitly with:

```bash
python3 tools/ingest_cantera_mechanism.py \
  --input mechanisms/h2o2_cantera.yaml --phase ohmech \
  --output mechanisms/h2o2_full.json \
  --module-name h2o2_full_mechanism_mod \
  --loader-name load_h2o2_full_mechanism \
  --kernel-name h2o2_full_production_rates \
  --jacobian-name h2o2_full_mass_fraction_jacobian \
  --thermo-loader-name load_h2o2_full_thermo_data \
  --transport-loader-name load_h2o2_full_transport_data \
  --symbol-prefix h2o2_full \
  --chemistry-integrator implicit \
  --description \
    "Full 10-species, 29-reaction H2/O2 bundle ingested from pinned Cantera YAML"
python3 tools/generate_elementary_mechanism.py \
  --input mechanisms/h2o2_full.json \
  --output src/generated/h2o2_full_mechanism_mod.F90
```

The importer publishes through a temporary sibling and atomic replacement;
it never exposes a partially validated JSON file. This command is a build-time
source update, not a runtime mechanism-selection interface.

To compile and execute a selected normalized bundle without modifying the
fixed application modules:

```bash
cmake -S . -B build-selected -DCMAKE_BUILD_TYPE=Release \
  -DPELEF_ENABLE_TESTS=OFF \
  -DPELEF_ENABLE_MPI=ON \
  -DCMAKE_Fortran_COMPILER=/path/to/gfortran \
  -DMPI_Fortran_COMPILER=/path/to/mpifort \
  -DMPIEXEC_EXECUTABLE=/path/to/mpiexec \
  -DPELEF_MECHANISM_BUNDLE="$PWD/mechanisms/h2o2_full.json" \
  -DPELEF_MECHANISM_SOURCE="$PWD/mechanisms/h2o2_cantera.yaml"
cmake --build build-selected \
  --target pelef_mechanism_probe pelef0d_selected \
    pelef_reactive_1d_selected pelef_amr_reactive_1d_selected \
    pelef_reactive_2d_selected pelef_reactive_eb_2d_selected \
    pelef_reactive_eb_amr_2d_selected pelef_reactive_3d_selected \
    pelef_reactive_eb_3d_selected \
    pelef_mpi_reactive_1d_selected \
    pelef_mpi_amr_reactive_1d_selected --parallel
./build-selected/pelef_mechanism_probe
./build-selected/pelef0d_selected cases/selected_reactor/h2o2_full.nml
./build-selected/pelef_reactive_1d_selected \
  cases/selected_reactive_1d/h2o2_full.nml
./build-selected/pelef_amr_reactive_1d_selected \
  cases/selected_amr_reactive_1d/h2o2_full_two_level_selected.nml
./build-selected/pelef_reactive_2d_selected \
  cases/selected_reactive_2d/h2o2_full_chemistry_plm.nml
./build-selected/pelef_reactive_eb_2d_selected \
  cases/selected_reactive_eb_2d/h2o2_full_transport_selected.nml
./build-selected/pelef_reactive_eb_amr_2d_selected \
  cases/selected_reactive_eb_amr_2d/h2o2_full_selected.nml
./build-selected/pelef_reactive_3d_selected \
  cases/selected_reactive_3d/h2o2_full_transport.nml
./build-selected/pelef_reactive_eb_3d_selected \
  cases/selected_reactive_eb_3d/h2o2_full_transport_chemistry.nml
/path/to/mpiexec -n 4 ./build-selected/pelef_mpi_reactive_1d_selected \
  mpi_reactive_1d_selected_np4.csv
/path/to/mpiexec -n 4 \
  ./build-selected/pelef_mpi_amr_reactive_1d_selected \
  cases/selected_mpi_sparse_amr_restart_1d/selected_reference.nml \
  selected_sparse_mpi_amr_np4.csv
```

The MPI launcher, wrapper, and `CMAKE_Fortran_COMPILER` must belong to one
toolchain. Remove foreign MPI include/library environment paths before
configuration; otherwise the wrapper can compile against one `mpi_f08` module
while the linker resolves `libmpi` from another implementation.

The optional stiff backend requires the pinned SUNDIALS `7.2.0` installation
to contain its official static Fortran-module libraries. C-only and
shared-only installs, and single- or extended-precision ABIs, are rejected.
After building that dependency with `BUILD_FORTRAN_MODULE_INTERFACE=ON`,
double precision, and static libraries enabled, use:

```bash
cmake -S . -B build-selected-cvode -DCMAKE_BUILD_TYPE=Release \
  -DPELEF_ENABLE_TESTS=OFF \
  -DPELEF_ENABLE_SUNDIALS=ON \
  -DCMAKE_PREFIX_PATH=/path/to/sundials-7.2.0-install \
  -DPELEF_MECHANISM_BUNDLE="$PWD/mechanisms/h2o2_full.json" \
  -DPELEF_MECHANISM_SOURCE="$PWD/mechanisms/h2o2_cantera.yaml"
cmake --build build-selected-cvode --parallel
./build-selected-cvode/pelef0d_selected \
  cases/selected_reactor/h2o2_full_cvode.nml
```

The baseline SUNDIALS commit and all recursive upstream revisions are frozen
in `references/pelec_baseline.json`. The CVODE build uses BDF with the serial
vector and dense matrix/linear-solver modules; it does not use a C-only
compatibility shim.

`PELEF_MECHANISM_SOURCE` is optional, but when supplied its basename and
SHA-256 must match the bundle provenance and a bundle selection is required.
Every complete schema-1 selected bundle must declare `chemistry_integrator`
as `explicit` or `implicit`; legacy kinetics-only generator inputs remain
exempt because they cannot configure a selected application.
The selected probe, `pelef0d_selected`, `pelef_reactive_1d_selected`,
`pelef_amr_reactive_1d_selected`, `pelef_reactive_2d_selected`,
`pelef_reactive_eb_2d_selected`, `pelef_reactive_eb_amr_2d_selected`, and
`pelef_reactive_3d_selected` are installed in a configuration that sets
`PELEF_MECHANISM_BUNDLE`.
`pelef_mpi_reactive_1d_selected` and
`pelef_mpi_amr_reactive_1d_selected`, and
`pelef_mpi_reactive_eb_patch_tree_2d_selected` are additionally installed
when `PELEF_ENABLE_MPI=ON`. The 0D reactor
namelist gives
`composition_count`, exact `composition_species` names, matching
`composition_mole_fractions`, time-integration controls, temperature,
pressure, and an output path. Bundle and source paths are CMake configure
dependencies, so changing either input and rebuilding refreshes all configured
applications and the source-hash expectation before regeneration.
`integrator` defaults to `implicit`; `cvode` is accepted only in a SUNDIALS
build. `maximum_steps` controls the native backend and
`cvode_max_internal_steps` controls CVODE's cumulative internal steps across
all requested output intervals. The selected regular-1D application reuses the
ordinary `reactive_1d` namelists, requires `chemistry_model = "selected"`, and
accepts the same exact-name composition arrays in `&reactive_1d`. The selected
AMR-1D application uses the same selected composition contract, requires
`amr_enabled = .true.`, and dispatches the established two-level, multilevel,
and dynamic-multipatch serial AMR modes. Checkpoint/restart remains outside
that selected front end. The selected
regular-2D application reuses `&reactive_2d`, requires
`chemistry_model = "selected"`, and supports the ordinary regular-grid 2D
periodic and physical-boundary controls. The selected single-level EB-2D
application additionally reuses `&embedded_boundary`, including plane/circle
geometry and embedded thermal/velocity wall controls, and rejects selected
AMR or checkpoint/restart use. It requires all four outer-domain boundaries
to be `outflow`; this is not general physical-boundary support. The selected
regular-3D application similarly reuses `&reactive_3d`, requires
`thermo_model = "selected"`, and remains periodic and single-level.
The selected regular-MPI-1D application deliberately accepts only an optional
output filename. It reuses the fixed 19-cell periodic verification setup,
bundle-order composition wave, native chemistry, molecular transport, and
ordered gather. Only the exact ordered H2/H two-species and pinned full-H2/O2
ten-species initialization profiles are admitted; it is not the general serial
`reactive_1d` namelist interface.
The selected sparse-MPI-AMR-1D application instead consumes the ordinary
`&reactive_1d` namelist, requires selected chemistry, AMR, and patch-tree mode,
and binds checkpoints to the bundle SHA-256, generated integrator, normalized
bundle-order composition, and original composite-integral baseline. A
checkpoint contains no MPI rank count or owner map, so a one-rank stop can be
continued on two or four ranks with deterministic redistribution. Use the
launcher selected by CMake/CTest when multiple MPI implementations are on
`PATH`.

With Ninja installed, presets cover ordinary, Cantera, and MPI builds:

```bash
cmake --preset debug
cmake --build --preset debug
ctest --preset debug
cmake --preset debug-cantera
cmake --preset mpi-release -DMPI_Fortran_COMPILER=/path/to/mpifort
```

## Representative runs

One-dimensional characteristic PLM:

```bash
./build/pelef cases/sod/sod_pelec_plm.nml
python3 tools/compare_sod.py --input sod_pelec_plm.csv
```

Two-dimensional periodic isentropic vortex:

```bash
./build/pelef2d cases/isentropic_vortex/vortex.nml
python3 tools/check_isentropic_vortex.py --input isentropic_vortex.csv
```

Three-dimensional periodic entropy wave:

```bash
./build/pelef3d cases/entropy_wave_3d/entropy_wave.nml
python3 tools/check_entropy_wave_3d.py \
  --input entropy_wave_3d.csv --nx 16 --ny 16 --nz 16 \
  --time 0.05 --maximum-l1 0.04
```

Passive two-species Sod problem:

```bash
./build/pelef_ms cases/multispec_sod/multispec_sod.nml
python3 tools/check_multispec_sod.py --input multispec_sod.csv
```

Elementary H2/O2 constant-volume reactor:

```bash
./build/pelef0d_h2o2 cases/zero_d_h2o2/reactor.nml
python3 tools/check_zero_d_h2o2.py --input zero_d_h2o2.csv
```

Reactive one-dimensional hotspot with characteristic PLM:

```bash
./build/pelef_reactive_1d cases/reactive_hotspot/hotspot.nml
python3 tools/check_reactive_hotspot.py --input reactive_hotspot.csv
```

The same case with monotone PPM and HLLC:

```bash
./build/pelef_reactive_1d cases/reactive_hotspot/hotspot_ppm.nml
python3 tools/check_reactive_hotspot.py --input reactive_hotspot_ppm.csv
```

The time-traced characteristic PPM path with the optional contact and shock
detectors enabled:

```bash
./build/pelef_reactive_1d \
  cases/reactive_hotspot/hotspot_characteristic_ppm.nml
python3 tools/check_reactive_hotspot.py \
  --input reactive_hotspot_characteristic_ppm.csv
```

Smooth general-EOS entropy wave:

```bash
./build/pelef_reactive_1d \
  cases/reactive_entropy_wave/entropy_wave.nml
```

General-EOS H2/N2 composition wave with HLLC:

```bash
./build/pelef_reactive_1d \
  cases/reactive_composition_wave/composition_wave.nml
python3 tools/check_reactive_composition_wave.py \
  --input reactive_composition_wave.csv
```


Dilute-gas transport coefficient probe:

```bash
./build/pelef_transport_probe transport_probe.csv
python3 tools/compare_transport_cantera.py --input transport_probe.csv
```

Periodic one-dimensional molecular-transport pulse:

```bash
./build/pelef_reactive_1d \
  cases/reactive_transport_1d/transport_pulse.nml
python3 tools/check_reactive_transport_1d.py \
  --input reactive_transport_pulse.csv --nx 96
```

Reactive two-dimensional hotspot with CTU and HLLC:

```bash
./build/pelef_reactive_2d cases/reactive_hotspot_2d/hotspot.nml
python3 tools/check_reactive_hotspot_2d.py \
  --input reactive_hotspot_2d.csv --nx 24 --ny 24
```

The same two-dimensional hotspot with characteristic PPM and the optional
contact/shock controls enabled:

```bash
./build/pelef_reactive_2d \
  cases/reactive_hotspot_2d/hotspot_characteristic_ppm.nml
python3 tools/check_reactive_hotspot_2d.py \
  --input reactive_hotspot_characteristic_ppm_2d.csv --nx 24 --ny 24
```

Oblique exact entropy-wave transport through the same 2D path:

```bash
./build/pelef_reactive_2d \
  cases/reactive_diagonal_wave_2d/diagonal_wave.nml
```

Oblique constant-pressure H2/N2 composition transport through the
characteristic-PPM/CTU path:

```bash
./build/pelef_reactive_2d \
  cases/reactive_diagonal_wave_2d/diagonal_composition_ppm.nml
```

Input-driven characteristic-PLM reactive EB hydro around a circular obstacle:

```bash
./build/pelef_reactive_eb_2d \
  cases/reactive_eb_circle_2d/uniform.nml
python3 tools/check_reactive_eb_circle_2d.py \
  --input reactive_eb_circle_2d.csv
```

Active-cell chemistry parity against the regular 2D path:

```bash
./build/pelef_reactive_eb_2d \
  cases/reactive_eb_chemistry_2d/reactive.nml
```

Runnable two-level reactive EB AMR hydrodynamics:

```bash
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_2d/uniform.nml
python3 tools/check_reactive_eb_amr_2d.py \
  --coarse reactive_eb_amr_coarse_2d.csv \
  --fine reactive_eb_amr_fine_2d.csv
```

Two-level reactive EB AMR thermal conduction and diffusive reflux:

```bash
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_transport_2d/reference.nml
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_transport_2d/transport.nml
python3 tools/check_reactive_eb_amr_transport_2d.py \
  --reference-coarse reactive_eb_amr_transport_reference_coarse.csv \
  --reference-fine reactive_eb_amr_transport_reference_fine.csv \
  --transport-coarse reactive_eb_amr_transport_coarse.csv \
  --transport-fine reactive_eb_amr_transport_fine.csv
```

Temperature-tagged conservative fine-patch movement:

```bash
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_2d/dynamic_hotspot.nml
python3 tools/check_reactive_eb_amr_dynamic_2d.py \
  --coarse reactive_eb_amr_dynamic_coarse_2d.csv \
  --fine reactive_eb_amr_dynamic_fine_2d.csv
```

Active-cell chemistry parity on both EB AMR levels:

```bash
./build/pelef_reactive_2d \
  cases/reactive_eb_amr_chemistry_2d/reference.nml
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_chemistry_2d/amr.nml
python3 tools/check_reactive_eb_amr_chemistry_2d.py \
  --reference reactive_eb_amr_chemistry_reference_2d.csv \
  --coarse reactive_eb_amr_chemistry_coarse_2d.csv \
  --fine reactive_eb_amr_chemistry_fine_2d.csv
```

Static three-level reactive EB AMR:

```bash
./build/pelef_reactive_2d \
  cases/reactive_eb_amr_chemistry_2d/reference.nml
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_three_level_2d/amr.nml
python3 tools/check_reactive_eb_amr_three_level_2d.py \
  --reference reactive_eb_amr_chemistry_reference_2d.csv \
  --root reactive_eb_amr_three_level_root_2d.csv \
  --middle reactive_eb_amr_three_level_middle_2d.csv \
  --finest reactive_eb_amr_three_level_finest_2d.csv
```

Static three-level checkpoint/restart parity:

```bash
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_three_level_restart_2d/reference.nml
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_three_level_restart_2d/checkpoint_stop.nml
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_three_level_restart_2d/restart.nml
python3 tools/check_reactive_eb_amr_three_level_restart_2d.py \
  --checkpoint three_level_restart.chk \
  --reference three_level_restart_reference_root.csv \
    three_level_restart_reference_middle.csv \
    three_level_restart_reference_finest.csv \
  --stopped three_level_restart_stopped_root.csv \
    three_level_restart_stopped_middle.csv \
    three_level_restart_stopped_finest.csv \
  --restarted three_level_restart_restarted_root.csv \
    three_level_restart_restarted_middle.csv \
    three_level_restart_restarted_finest.csv
```

Tag-driven dynamic parent and finest patches:

```bash
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_three_level_dynamic_2d/hotspot.nml
python3 tools/check_reactive_eb_amr_three_level_dynamic_2d.py \
  --root three_level_dynamic_root.csv \
  --middle three_level_dynamic_middle.csv \
  --finest three_level_dynamic_finest.csv
```

Reacting fine-to-root checkpoint/restart parity:

```bash
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_restart_2d/reference.nml
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_restart_2d/checkpoint_stop.nml
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_restart_2d/restart.nml
python3 tools/check_reactive_eb_amr_restart_2d.py \
  --checkpoint reactive_eb_amr_restart.chk \
  --reference reactive_eb_amr_restart_reference_coarse_2d.csv \
  --stopped reactive_eb_amr_restart_stopped_coarse_2d.csv \
  --restarted reactive_eb_amr_restarted_coarse_2d.csv \
  --fine reactive_eb_amr_restart_reference_fine_2d.csv \
    reactive_eb_amr_restart_stopped_fine_2d.csv \
    reactive_eb_amr_restarted_fine_2d.csv
```

Reacting two-child checkpoint/restart parity:

```bash
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_multipatch_restart_2d/reference.nml
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_multipatch_restart_2d/checkpoint_stop.nml
./build/pelef_reactive_eb_amr_2d \
  cases/reactive_eb_amr_multipatch_restart_2d/restart.nml
ctest --test-dir build --output-on-failure \
  -R '^regression_reactive_eb_amr_multipatch_restart_2d_'
```

## Project records

- [Completion roadmap](docs/completion_roadmap.md)
- [Porting plan](docs/porting_plan.md)
- [Architecture](docs/architecture.md)
- [Numerical methods](docs/numerical_methods.md)
- [PeleC responsibility mapping](docs/pelec_mapping.md)
- [State variables](docs/state_variables.md)
- [Parity strategy](docs/parity_strategy.md)
- [Implementation status](docs/implementation_status.md)
- [Design decisions](docs/design_decisions/)


### PeleF 0.17.0 two-dimensional transport example

```bash
./build/pelef_reactive_2d cases/reactive_transport_2d/transport_hotspot.nml
python3 tools/check_reactive_hotspot_2d.py \
  --input reactive_transport_hotspot_2d.csv --nx 20 --ny 20
```


### PeleF 0.18.0 physical-boundary examples

```bash
./build/pelef_reactive_2d cases/reactive_boundaries_2d/couette.nml
./build/pelef_reactive_2d cases/reactive_boundaries_2d/thermal_channel.nml
./build/pelef_reactive_2d cases/reactive_boundaries_2d/inflow_outflow.nml
./build/pelef_reactive_2d \
  cases/reactive_boundaries_2d/prescribed_species_wall.nml
```

Solid walls are species-impermeable unless `wall_species_* = "prescribed"`
supplies a zero-sum species mass-flux vector in kg/(m2 s), positive from the
wall into the gas. Slip walls remove tangential viscous stress; no-slip walls
reflect velocity about a prescribed wall velocity. Prescribed fluxes provide
a catalytic-wall transport interface, not a surface-reaction-rate model.


### Full pressure-dependent H2/O2 chemistry

Set `chemistry_model = "full_h2o2"` to use the 10-species, 29-reaction third-body/Troe mechanism with the implicit cell reactor.


### PeleF 0.24.0 MPI 1D verification

Configure with `-DPELEF_ENABLE_MPI=ON`, then run the MPI verification drivers
with 1, 2, 4, or 8 ranks. The foundation and multispecies drivers use 257 cells
so the block decomposition is intentionally uneven; the chemistry, transport,
and coupled-reactive drivers use smaller non-divisible workloads to exercise the
same decomposition and ordered-gather logic.
