# Distributed nodal EB gas/spray example

Build with `PELEF_ENABLE_MPI=ON` and `PELEF_ENABLE_SPRAY=ON`, then run from
this directory (choose a fresh directory for each run):

```sh
mpiexec -n 2 /path/to/pelef_mpi_eb_spray_3d tree.nml spray.nml
```

The fixed executable uses the existing elementary H2/O2 gas initialization.
A bundle configured with an MPI EB tree executable also builds
`<selected_mpi_eb_tree_target>_spray`. It accepts the same two paths and selected
composition fields in `/eb_tree/` (for example pure N2 via `composition_count=1`,
`composition_species(1)='N2'`, `composition_mole_fractions(1)=1.`).

This fixture has a 4 mm cube, a ratio-two left child, a planar nodal EB at
y=1.25 mm and two initially injected water parcels directed toward the wall.
The 8³ finest nodal file is required, including its 9³ nodes in x-fastest order.
Gas chemistry, all molecular transport terms, hydro and conservative particle
feedback advance together for 20 ns. Gas controls use the existing `/eb_tree/`
namelist. `/eb_spray/` accepts `nc`, `number`, component property arrays
(`vapor_names`, `densities`, `heat_capacities`, `reference_temperatures`,
`reference_pressures`, `latent_heats`, `fractions`), `total_mass`, `center`,
`axis`, `speed`, `half_angle`, `diameter`, `liquid_temperature`, and
`particle_options` derived-type members. Up to 16 liquid components and 125000
initial records are allowed, further limited by `particle_options%maximum_parcels`.
`injection_events` (1--1024), `injection_times` and `injection_scales` select
ordered events; defaults are one event at zero with unit mass scale. Scales
multiply parcel multiplicity, preserving drop size and composition. Equal times
are allowed; zero scales consume an event but no IDs. The time step stops at
each event, including final-time injection; later events remain pending. A
nondecreasing finite schedule and the global live-population cap are enforced.
For example, `injection_events=3,injection_times=0.,1.e-8,2.e-8,
injection_scales=1.,.5,2.` injects at the start, midpoint and end.

Output:

- `result.csv`: exposed composite gas cells, using the existing MPI EB writer.
- `result.csv.budget.csv`: combined gas-plus-particle initial/final inventories
  and exterior, wall, chemistry and injected exchanges, all with change-in-system signs.
  Initial inventory is gas before injection, including for a zero-time event.
- `result.csv.particles.rank<R>.csv`: local parcel IDs, position, velocity,
  mass, temperature, multiplicity and per-drop component masses. Empty ranks
  write a header. IDs are global; concatenate and sort by ID for comparison.
  Existing particle CSV files are protected and cause the run to fail.

This first executable connection supports a fixed hierarchy and scheduled cone
injection. Regridding and child creation are explicitly rejected. It inherits root-contained specular EB drift and the
composite stencil approximation. It is a synthetic integration fixture, not
nozzle, diesel ignition, convergence or physical validation.

## Checkpoint and continuation

Set `checkpoint_file='saved'`, `checkpoint_interval_steps=2` and optionally
`stop_after_checkpoint=.true.` in `/eb_tree/`. Each accepted checkpoint step
publishes an immutable `saved.step<N>` manifest, plus `.gas` and `.particles`
payloads. `.pending` files reserve the names and remain beside successful
archives. Existing files are never replaced. Retain the complete set together.

Resume with `restart_file='/absolute/path/saved.step2'`, a fresh output path,
and the same tree/spray model inputs. The MPI rank count may change. Final time,
step limit and output/checkpoint controls may change within the existing gas
restart contract. A checkpoint is taken after the accepted step and before
injections due at its new time; resuming processes those events exactly once.
All particle fields, event position, ID watermark, injected inventory, gas state
and cumulative coupled ledgers continue together.

The context binds normalized spray namelist values and original nodal values;
the existing gas archive validates its mechanism/configuration and geometry.
A mismatched context, incomplete payload or inconsistent metadata rejects the
restart before publishing the staged state. The manifest appears only after
both payloads have closed successfully. This provides atomic publication on the
shared filesystem, not power-loss durability or a cross-compiler format guarantee.
The serialized spray/nodal context is capped at 130000 characters; larger inputs
are rejected when checkpointing is requested.

## Whole-patch load redistribution

`rebalance_interval_steps=N` samples the existing gas-kernel work counters every
N accepted steps. `rebalance_minimum_gain` is the required predicted gain
(default 0.1); zero accepts any positive gain. Gas/temperature and particles with
their component masses move atomically to the proposed patch owners. If particle
migration fails, the original ownership and all state remain intact. A move is
reported as `Rebalanced at step:`; an already balanced layout may stay unchanged.
This uses measured gas cost, not a separate model of particle work. The
partition is not part of the restart identity, so continuation can redistribute
again under a different rank count. The geometry and hierarchy do not change.
