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
PCM/characteristic-PLM face reconstruction, physical boundaries, chemistry and
molecular transport. Existing validation requires at least four root cells in
each direction. The application is also registered for installation.

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

Rank zero gathers state once for the initial inventory and gathers state/temperature
at final output and requested checkpoints. Final CSV rows and budget format match
serial. Numerical time stepping does not gather full fields.
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
format as serial. Rank zero reads/writes the checkpoint, validates its context,
scheduled/adaptive geometry, consumed event or adaptation-check count, EOS, averages and ledger, then redistributes payloads according
to the current rank count. Fixed/selected checkpoints can be exchanged with the
corresponding serial frontend for static, scheduled or automatic runs in either timestep mode. The selected bundle
SHA-256 remains part of checkpoint identity. A requested checkpoint stop returns
before final CSV output; a restart resumes step counters without resaving the
loaded step. Checkpoint files must not alias input, schedule, any scheduled topology or output paths. Input
and output lexical aliases are rejected after resolving relative paths against
the working directory. This POSIX path check does not resolve symlinks/hard links.
Whole-patch ownership and root output are not a production scaling claim.

## Remaining scope

This frontend supports static, scheduled or temperature-adapted topologies with serial ownership.
The MPI frontend supports static/scheduled/automatic execution and restart with
common or recursive timesteps. Characteristic-PLM hydro enables limited
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
The same schedule/checkpoint controls, output alias checks and root CSV gathers as
the fixed MPI frontend apply, including recursive subcycling.

Automatic MPI restart retains the saved temperature guesses, tile geometry and
adaptation-check counter. Rank zero validates the configured depth, ratios, tile
alignment and node budget before broadcasting geometry metadata and scattering
payloads to the current owners. The saved step is not tagged again on restart.
Changing adaptation thresholds or cadence changes checkpoint identity and is rejected.
