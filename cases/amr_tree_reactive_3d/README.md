# Experimental serial regular 3D patch-tree application

Build and run from the repository root:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_amr_tree_reactive_3d
mkdir -p build/tree-example
cd build/tree-example
../spray-mpi-migration-gnu-debug/pelef_amr_tree_reactive_3d \
  ../../cases/amr_tree_reactive_3d/hotspot.nml \
  ../../cases/amr_tree_reactive_3d/branched.dat
```

The executable accepts `<reactive_3d.nml> <topology.dat> [schedule.dat]`. Flow input uses the
existing `reactive_3d` configuration, including elementary/full-H2/O2 models,
PCM/characteristic-PLM/characteristic-PPM face reconstruction, physical boundaries, chemistry and
molecular transport. Existing validation requires at least four root cells in
each direction. The application is also registered for installation.

## SDC

With `PELEF_ENABLE_SPRAY=ON`, the serial/MPI fixed and selected frontends also route
LES into SDC. Add this group to the flow input:

```fortran
&amr_tree_les
  enabled=.true., model='smagorinsky', cs=.17, prt=.7, sct=.7, ci=0
/
```

Models are `smagorinsky`, `wale`, `vreman` and `dynamic`. Dynamic-only controls
are `test_model='box'`, `test_ratio=2`, `grid_ratio=1`, and
`harmonic_mean=.false.`. LES requires `time_integrator='sdc'` and supports the
existing common/recursive setting. A conservative diffusion limit of 0.1 uses
the current finest cell width, including in recursive mode. LES settings are
part of checkpoint identity; changing or removing enabled LES rejects restart.
Omitting the group leaves LES disabled. The selected route is wired but the
[focused application check](../../docs/validation/les-mpi-amr-app-20260927.md)
was run only for the fixed mechanism. Serial/MPI output and checkpoint
interchange also have [focused evidence](../../docs/validation/les-serial-amr-app-20260927.md).
Particle-SDC routing remains separate work.


Serial and MPI fixed/selected frontends accept `time_integrator='sdc'` with
`characteristic_ppm` and the `pelec` Riemann solver. Run `sdc.nml` with `sdc.dat`
for a small two-level reacting example with transport and external forcing.
Reaction-source estimates persist across successful steps. The SDC budget CSV
adds `source_change`, separate from `boundary_change`.

This route uses a common timestep by default. Add
`&amr_tree_timestepping subcycling=.true. /` in serial or MPI runs to advance finer levels with
refinement-ratio substeps and synchronize their accumulated fluxes. SDC checkpoint/restart preserves reaction
estimates and external-source inventories in V3 (plain) or V4 (indexed) records.
Scheduled and automatic regridding transfer signed reaction estimates by the
same volume-overlap remap used for the gas. All three fields (gas, temperature,
reaction estimates) publish together after successful remapping. Recursive SDC uses the same checkpoint fields and binds the subcycling choice
to restart context. The MPI checkpoint library reads/writes these SDC formats,
including changed-rank source redistribution. MPI applications support common-step and recursive
SDC with checkpoint/restart and regridding. Rebalancing packs gas and reaction
estimates together for owner migration. Recursive MPI SDC retains owner-local state and accumulates child fluxes before
coarse/fine synchronization.
Split stepping retains its existing capabilities and CSV format; nonzero
external forcing requires SDC.

## Pressure outlets

The serial and MPI fixed/selected tree frontends also accept `inflow_mach(6)` in
`reactive_3d`. Negative values retain reservoir inflow. A value from 0 to 100
activates Mach inflow on that `inflow` face, with positive
`inflow_temperature` and zero `inflow_velocity`. `inflow_pressure=0` extrapolates
interior pressure; a positive value prescribes pressure in Pa. Composition comes
from the current interior donor. Each step limits hydro and molecular transport
using the generated ghosts, with the existing refinement scaling for subcycling.
Checkpoint context binds these controls. MPI computes the minimum bound over
owned nodes and reduces it across ranks, including ranks without any nodes.

The serial and MPI fixed/selected frontends accept `outlet_pressure(6)` in `reactive_3d`.
For example, `outlet_pressure(2)=102325` sets a 102325 Pa target on an `outflow`
x-high face; zero retains extrapolation. Targets must be finite, positive only
on outflow faces, and at most 1e15 Pa. The pressure ghost extrapolation feeds
both hydro and molecular transport with common or recursive timesteps.
Nonphysical ghost states reject the candidate step and trigger the existing
bounded timestep retry. This does not specify a backflow reservoir or prove
nonreflection.

Targets are part of checkpoint context: changing or removing them at restart
is rejected. Scheduled/automatic remapping preserves volume inventories and
does not sample physical ghost states. MPI stepping checks target agreement
collectively. Rank-changing and serial/MPI restart retain the same pressure
identity requirements.

## Configure-time selected mechanisms

When `PELEF_MECHANISM_BUNDLE` is configured, the bundle generator also provides
`pelef_amr_tree_reactive_3d_selected`. It accepts the same topology and optional
schedule arguments. Set `thermo_model='selected'` and supply species names and
mole fractions through `composition_count`, `composition_species`, and
`composition_mole_fractions`, matching the selected bundle.

The selected frontend loads generated thermo/reactions/transport, validates the
initial/inflow/isothermal-wall temperature envelope and forwards the generated
chemistry integrator. It logs the bundle SHA-256 and integrator. Its dedicated
runtime links selected libraries without `pelef_core` or fixed database imports.
For custom bundle targets, `pelef_add_mechanism_bundle` accepts
`AMR_TREE_REACTIVE_3D_NAME`. The selected frontend is included in the bundle's
normal install target set; an installed-binary audit is separate evidence.

## Topology file

The first record is the number of child patches (zero is a root-only run).
Each following record contains:

```text
parent ilower iupper jlower jupper klower kupper ratio
```

Root node is 1. Record `i` defines node `i+1`; its parent must already exist.
Bounds are inclusive parent-cell indices. Siblings cannot overlap. Refinement
ratios must agree at each depth. Counts/depths are dynamically allocated, with
checked coordinate/extent limits. `branched.dat` defines two adjacent root
children and one grandchild: four nodes over three levels.

Each node samples the analytic initial condition at its own global cell centers.
Covered parents are then restricted and temperatures recovered. Time integration
uses synchronous R/2-T/2-H-T/2-R/2, with disabled operators skipped. The timestep
includes all-node hydro/transport bounds, touched isothermal walls and inflow
states at the finest spacing. Rejected steps halve the timestep up to 20
attempts; only successful steps advance time and boundary budgets.

## Scheduled regridding

An optional third argument gives a schedule. The first record is its event count;
each following record is `time 'topology-path'`. Times must increase strictly,
be positive, and not exceed `final_time`. Quote paths, including absolute paths;
relative paths resolve beside the schedule file. All topology inputs are read
before evolution. `regrid.dat` moves/refines the patches with ratio 3, removes
all children, then restores `branched.dat` at final time:

```text
3
2.5e-10 'moved_ratio3.dat'
5.0e-10 'root.dat'
1.0e-9 'branched.dat'
```

Time steps stop exactly at each event. Conserved overlap projection and EOS
recovery occur before the new topology is published. An absolute-cell-inventory
scaled conservation check must pass at `5e-12`; boundary budgets and original
initial inventories are carried through unchanged. Final-time events run before
output. Reservoir bounds use the finest spacing across the initial and scheduled
trees; per-step hydro/diffusion/wall bounds use the active tree. Logs report event
times, node/level counts, and maximum remap error. This is a prescribed schedule,
not automatic solution-based tagging.

## Hierarchy subcycling

The serial and MPI fixed/selected frontends accept:

```fortran
&amr_tree_timestepping
  subcycling=.true.
/
```

The default is `.false.` (the original common timestep). With subcycling enabled,
one accepted application step is a coarse R/2-T/2-H-T/2-R/2 interval. Hydro uses
one child step per refinement ratio, while each transport half uses the ratio
squared. Chemistry operates on exposed cells over the two coarse half-intervals.
Every completed fine interval returns corrected fluxes to its parent. The entire
coarse interval rolls back if any operator fails; the existing timestep-halving
retry then applies.

Timestep selection uses scaled hydro, diffusion and isothermal-wall bounds.
Reservoir bounds use root spacing because the corresponding spatial ratios are
cancelled by the hydro/transport substep counts. Transport bounds conservatively
limit the full coarse dt even though each transport solve spans dt/2. Regrid
cadence, `maximum_steps` and checkpoint intervals count coarse steps. Prescribed
event times still clip the coarse step exactly; regridding and checkpoint writes
occur only after the full interval succeeds. Static, prescribed and automatic
geometries use the same dispatch. Checkpoint context includes subcycling mode;
changing it at restart is rejected. No checkpoint contains an unfinished substep.

## Automatic temperature adaptation

Instead of a schedule argument, add this optional group to the input file:

```fortran
&amr_tree_adaptation
  regrid_interval_steps=1
  refinement_levels=2
  refinement_ratios=2,2
  temperature_thresholds=1150.0,1200.0
  tile_shape=2,2,2
  tag_buffer=1
  maximum_nodes=4096
/
```

The initial topology still comes from the second command-line argument. After
each requested number of successful steps, exposed-cell temperatures select
parent-cell tiles for refinement; conservative remapping applies any geometry
change. This includes a final-time update. The initial state is not retagged at
step zero. A check that selects unchanged geometry increments the regrid counter
without remapping. The logged `Regrids` count is therefore the number of checks.

`refinement_levels` permits 0–16 levels above root (zero coarsens to root).
Positive thresholds must not decrease with depth; ratios must be at least 2.
Tiles clip at parent edges; the nonnegative buffer is measured in parent cells
at each depth and wraps periodic directions. `maximum_nodes` includes root;
exceeding it rejects the proposal instead of silently dropping tags. Initial
proposal validation occurs before stepping, and reservoir timestep bounds cover
the finest permitted spacing. Omit the group or set `regrid_interval_steps=0`
to disable adaptation. Prescribed schedules and automatic adaptation are separate
modes and cannot be combined.

Fixed and selected mechanisms share these controls. Automatic restart preserves
the last generated topology; it does not retag saved temperatures, since remapping
may have changed the temperatures that originally triggered refinement. Controls
are part of checkpoint identity. The reader checks the consumed update count and
permitted depth, ratios, tile alignment/extents and node budget. These checks do
not reconstruct historical tag decisions. New updates begin at the next due
successful step. Both optional namelist group headers must start new lines; their
order does not matter.

## Checkpoint and continuation

Both fixed and selected frontends accept an optional group in the input file
(group header starts a new line):

```fortran
&amr_tree_checkpoint
  checkpoint_file='tree.chk'
  checkpoint_interval_steps=1
  stop_after_checkpoint=.true.
/
```

The interval counts successful steps. A checkpoint is written after any regrid
at that time, including a final-time regrid. Stopping returns successfully before
writing final CSVs; existing CSVs remain untouched. Omit the group to disable
checkpointing. To continue, use the same physical input, initial topology and
full schedule, replacing the group with:

```fortran
&amr_tree_checkpoint
  restart_file='tree.chk'
/
```

Paths resolve from the working directory. Final time, maximum steps and output
path may change; all physics, mechanism records/identity, initial topology and
ordered schedule contents must match. Topology filenames may change if their
contents do not. The reader checks consumed events, active geometry, run limits,
initial inventories and the conserved boundary ledger before continuing. Saved
temperatures and accumulated diagnostics are restored exactly. Restart never
replays a consumed event. Saving can also remain enabled during continuation;
the restored step is not written again. Checkpoint and restart paths may be the
same, but must not alias configuration/topology/schedule or CSV paths. The path
guard is lexical; symlink and hard-link aliases are not resolved.

The snapshot writer validates data before opening the destination. An I/O failure
can still leave a partial checkpoint. Exact context text is compiler/runtime
sensitive; cross-compiler portability has not been qualified.

## Output

`output_file` is a single CSV of **exposed cells only**, so covered parent cells
are not counted twice. `node,level,dx,dy,dz` precede the standard reactive-3D
columns. Cell volume is `dx*dy*dz`; field integrals must use that weight.
The companion `<output_file>.budget.csv` records `component,initial,final,
boundary_change`. Components are density, x/y/z momentum, total energy, then
species densities in mechanism order. Chemistry changes individual species;
the corresponding source contribution is not a boundary flux.

`hotspot.nml` exercises elementary chemistry, transport and moving no-slip walls.
The small case produces 127 exposed cells and is an implementation regression,
not physical validation. Output write failure may leave a partial CSV; the
files are not checkpoints.

## Static MPI execution

With an MPI-enabled build:

```sh
cmake --build build/tree-distribution-mpi-debug --target pelef_mpi_amr_tree_reactive_3d
mpiexec -n 3 build/tree-distribution-mpi-debug/pelef_mpi_amr_tree_reactive_3d \
  cases/amr_tree_reactive_3d/hotspot.nml cases/amr_tree_reactive_3d/branched.dat
```

This fixed-mechanism executable uses the same `reactive_3d` input and topology
format, supports `elementary` and `full_h2o2` loading, and assigns whole nodes
round-robin by node index. Each rank initializes and advances only owned payloads;
empty ranks are supported. It uses common timesteps and the full reaction/transport/
hydro composition with up to 20 timestep-halving attempts.

Initial and final inventories use serial-order compensated small-vector transfers.
Restart ledger validation obtains its initial absolute-inventory scale the same
way, without gathering initialized fields. Final CSV rows are written
by patch owners through MPI-IO at predetermined byte offsets, preserving serial
node/cell order and formatting. All rows are validated before opening the file;
owner buffers hold bounded chunks. The final budget uses serial-order compensated
accumulators with only small-vector broadcasts, and rank zero writes its CSV.
Initialization, numerical stepping, final CSV output and checkpoint saves do not
gather full fields. Checkpoint saves use owner-local MPI-IO after distributed
EOS/parent-average validation; the bytes remain compatible with serial. Checkpoint
loading indexes V1 cell records on rank zero with one-cell scratch space, then
each owner reads only its assigned node ranges. The index scan remains serial;
no complete root state is constructed. I/O errors after opening may leave a
partial output; this is not an atomic filesystem transaction.
The bounded regression checks byte equality for the elementary hotspot case.

An optional third CLI argument supplies the same regrid schedule as serial.
Event times clip the timestep exactly. Each event conservatively remaps old
exposed cells to the new round-robin owners, checks component inventories and
updates geometry before checkpointing/output. The inventory check transfers
compensated accumulators in serial cell order without gathering field payloads.
The `amr_tree_adaptation` group accepts the same controls as serial and is exclusive
with a schedule argument. After each configured step interval, owners identify
hot exposed cells and generate the same tile geometry on every rank. Geometry
changes use sparse conservative remapping before checkpoint/output. The requested
finest spacing bounds common-step reservoir timesteps; recursive mode uses root
spacing. The `amr_tree_timestepping` group selects the recursive mode described above. The `amr_tree_checkpoint` group uses the same controls and
format as serial. Owners write their checkpoint rows directly. On restart, all
ranks read the shared header; rank zero streams through the records to build a
node-position index and check the footer. Owners seek to their assigned ranges.
Context, scheduled/adaptive geometry, consumed event or adaptation-check count,
EOS, averages and ledger are validated before replacing the current fields. Fixed/selected checkpoints can be exchanged with the
corresponding serial frontend for static, scheduled or automatic runs in either timestep mode. The selected bundle
SHA-256 remains part of checkpoint identity. A requested checkpoint stop returns
before final CSV output; a restart resumes step counters without resaving the
loaded step. Checkpoint files must not alias input, schedule, any scheduled topology or output paths. Input
and output lexical aliases are rejected after resolving relative paths against
the working directory. This POSIX path check does not resolve symlinks/hard links.
Whole-patch ownership and distributed CSV output are not a production scaling claim.

### Indexed checkpoints

Add `checkpoint_indexed=.true.` to `amr_tree_checkpoint` when saving to select
V2 with checked relative byte offsets for every node and the payload end. Serial
and MPI fixed/selected frontends support this option and produce identical bytes.
Readers detect V1/V2 automatically; the format choice does not affect physical
checkpoint identity. The default remains V1.

V2 lets MPI readers seek directly to owned ranges and the footer without scanning
all cell records on rank zero. Its canonical fixed-width ASCII records use LF
newlines; the stored index is checked against geometry, variable count and actual
consumed record lengths. V1 retains the flexible list-directed text path and its
bounded-memory index scan. Older executables that only recognize V1 cannot read V2.

## MPI work-weighted patch ownership

The MPI frontends accept this optional group:

```fortran
&amr_tree_distribution
  partitioning='weighted'
/
```

The default is `round_robin`. Weighted assignment places the largest estimated
patch workload on the least-loaded rank, breaking ties by node/rank index.
Weights are cell counts for common stepping, cell counts times cumulative
refinement ratios for recursive hydro, and times squared ratios when recursive
transport is enabled. These are cell-step estimates, not measured chemistry or
communication costs. Whole patches remain indivisible.

The policy is applied at initialization, regridding and checkpoint reload.
It is not part of physical checkpoint identity: rank count and assignment policy
may change at restart, including interchange with serial. Regridding transfers
payloads directly to the new owners through the existing sparse remap.

Optional runtime rebalancing is configured in the same group:

```fortran
&amr_tree_distribution
  partitioning='weighted'
  rebalance_interval_steps=10
  rebalance_minimum_gain=0.1
/
```

The interval defaults to zero (disabled). On each positive multiple of the interval,
the application measures that completed coarse step's owner-local hydro/transport
flux and update work plus chemistry solves, including recursive substeps. Failed
attempts are discarded. It then migrates U/T only if the predicted maximum rank
cost falls by at least the minimum fractional gain (default 0.1, valid range [0,1)).
Geometry is unchanged and values transfer exactly. `Rebalance step` logs each
decision. Timings are fresh after regrid/restart, are not checkpointed, and do not
enter physical restart identity; cadence and gain can change on continuation.

Wall-clock noise can change owner choices. Costs omit communication, restriction,
EOS recovery outside chemistry, timestep selection and migration overhead; this
is not a demonstrated speedup. Whole-patch constraints still apply.

## Remaining scope

This frontend supports static, scheduled or temperature-adapted topologies with serial ownership.
The MPI frontend supports static/scheduled/automatic execution and restart with
common or recursive timesteps. Characteristic-PLM and characteristic-PPM hydro enable limited
linear coarse donors with the configured limiter for both hydro and transport.
PCM reconstruction retains PCM donors. Common timesteps and recursive subcycling are available; globally second-order
AMR convergence and large-tree scaling remain unqualified.

### Selected-mechanism MPI application

With `PELEF_ENABLE_MPI=ON` and a schema-1 `PELEF_MECHANISM_BUNDLE`, build
`pelef_mpi_amr_tree_reactive_3d_selected` and run:

```sh
mpiexec -n 3 ./build/selected/pelef_mpi_amr_tree_reactive_3d_selected input.nml topology.dat
```

Use `thermo_model='selected'` and the bundle's named mole-fraction composition.
The frontend validates the mechanism and initial temperature range, logs its
SHA-256/integrator on rank zero, and forwards that integrator to owner-local
stepping. Its isolated MPI runtime links selected libraries without `pelef_core`
or fixed mechanism databases. Custom bundle targets can use
`MPI_AMR_TREE_REACTIVE_3D_NAME`; the target participates in bundle installation.
The same schedule/checkpoint controls, output alias checks and distributed CSV writes as
the fixed MPI frontend apply, including recursive subcycling.

Automatic MPI restart retains the saved temperature guesses, tile geometry and
adaptation-check counter. Rank zero validates the configured depth, ratios, tile
alignment and node budget before broadcasting geometry metadata and scattering
payloads to the current owners. The saved step is not tagged again on restart.
Changing adaptation thresholds or cadence changes checkpoint identity and is rejected.

## Characteristic PPM

Set `reconstruction='characteristic_ppm'` in `&reactive_3d` for fixed or
selected serial/MPI tree gas applications. PPM uses four ghost layers and
one time-centered CTU step per level interval, with flattening enabled and
contact steepening disabled. Both common and recursive timesteps use
integrated PPM fluxes for reflux. Molecular transport uses limited coarse
donors, and Mach inlet timestep evaluation uses four-layer donors.

Checkpoint context includes the reconstruction selector; switching to PCM
or PLM on continuation is rejected. The
[focused application checks](../../docs/validation/ppm-amr-application-20260925.md)
cover reacting molecular transport, serial/MPI output equality and fixed
mechanism rank-changing continuation. This does not establish convergence
order or exact PeleC primitive-pressure parity.

### Artificial viscosity

`difmag=-1` (default) selects 0.1 with PPM/pelec and zero otherwise. Explicit
nonnegative values work with PCM, PLM and PPM; zero disables
the correction. The resolved coefficient reaches every gas hydro level and is
stored in checkpoint context. Automatic and equivalent explicit settings can
resume each other. Changing the effective coefficient rejects restart. Older
PPM checkpoints that stored the unresolved -1 sentinel are incompatible; default
PCM/PLM context is unchanged, including equivalent explicit zero. See
[coefficient identity checks](../../docs/validation/nonppm-difmag-input-20260925.md)
and [PPM application restart evidence](../../docs/validation/difmag-amr-input-20260925.md).

## SDC reactor selection

With `PELEF_ENABLE_SUNDIALS=ON`, fixed and selected serial/MPI tree applications
accept `sdc_reactor='cvode'` in `&reactive_3d`. It requires
`time_integrator='sdc'` and `chemistry_enabled=.true.`. The default `native`
keeps the existing forced reactor. CVODE is used for each cell's forced
reaction solve in both common and recursive stepping, including the LES
application wrapper. This selector is independent of the selected bundle's
split-step chemistry integrator.

The reactor selector is part of the existing full-configuration checkpoint
identity. Continue with the same selector; changing it rejects the checkpoint.
Builds without SUNDIALS reject `cvode` during configuration input.
