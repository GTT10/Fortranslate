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
- `result.csv.budget.csv`: combined gas, particle and film initial/final inventories
  and exterior, wall, chemistry and injected exchanges, all with change-in-system signs.
  Initial inventory is gas before injection, including for a zero-time event.
- `result.csv.particles.rank<R>.csv`: local parcel IDs, `material_id` (injection port), position, velocity,
  mass, temperature, multiplicity and per-drop component masses. Empty ranks
  write a header. IDs are global; concatenate and sort by ID for comparison.
  Existing particle CSV files are protected and cause the run to fail.

This executable supports scheduled cone injection and the existing EB tree
regrid controls. KH child shedding remains explicitly rejected; optional wall impacts create splash children. It inherits root-contained EB flight and the
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
All particle and film fields, event position, ID watermark, injected inventory, gas state
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
(default 0.1); zero accepts any positive gain. Gas/temperature, particles and films with
their component masses move atomically to the proposed patch owners. If particle or film
migration fails, the original ownership and all state remain intact. A move is
reported as `Rebalanced at step:`; an already balanced layout may stay unchanged.
This uses measured gas cost, not a separate model of particle work. The
partition is not part of the restart identity, so continuation can redistribute
again under a different rank count. Rebalancing changes ownership only; the physical geometry remains fixed.

## Regridding and restart stages

The existing `/eb_tree/` controls `regrid_count`, `regrid_times`,
`regrid_child_count`, `regrid_parent`, `regrid_ratio`, `regrid_lower` and
`regrid_upper` now remap gas, particles and films together. Refinement may move,
coarsen to the root or grow again on the shared refinement lattice. Supply the
original nodal field at the maximum lattice resolution required by the complete
schedule. Coarsening retains that field for collision/visibility queries.
The existing `temperature_adaptive` tagger also routes through the coupled remap;
its thresholds, depth and capacity controls are unchanged.

At a step ending on a regrid time, the hierarchy changes before balance checks
and checkpointing. Injection events due at that time run at the beginning of
the next iteration. Restart validates the saved regrid schedule/stage and builds
a private target hierarchy before restoring gas and particles. The remaining
regrid and injection events then continue. The focused application evidence
covers prescribed schedules; adaptive tagging is not newly qualified by it.

## Independent injection ports

`port_count` selects 1--16 ports; it defaults to one. `event_ports(e)` selects
the port for each entry in the existing globally ordered injection schedule
(default one). Equal times may select different ports; each event's scale
multiplies that port's base injected mass. Ports may have different particle
counts. Positive events reserve only their own count of consecutive global IDs;
zero-scale events reserve none.

Per-port arrays are `port_numbers`, `port_centers(3,16)`, `port_axes(3,16)`,
`port_masses`, `port_speeds`, `port_angles`, `port_diameters`,
`port_temperatures`, and `port_fractions(16,16)` (component, port). Unspecified
entries inherit the corresponding legacy single-port input. The internal unset
marker is -1 for scalar/fraction arrays and `huge(real64)` for complete vector
columns. A supplied fraction column must specify all active fractions, zero
unused entries, and sum to one. All ports share the component property table;
their compositions, cone geometry and injection state can differ.

For example, add a second port and select it at the second event:

```fortran
port_count=2,event_ports=1,2,
port_numbers(2)=3,port_centers(:,2)=.003,.001250001,.0015,
port_masses(2)=2.e-10,port_speeds(2)=30.,
port_diameters(2)=2.e-5,port_temperatures(2)=360.,
```

The full normalized port configuration and event selections are bound into the
checkpoint context. CSV `material_id` records the port number. The old single
port input remains valid; checkpoints produced before this input/context extension
are rejected by the strict context comparison. Events are specified in one
ordered table; separate per-port rate-curve files are not connected here.

## Port mass-flow tables

Set `port_profile_points(p)` to 2--64 and provide
`port_profile_times(:,p)` (seconds) and `port_profile_rates(:,p)` (kg/s).
Times must strictly increase; rates are finite and within 0--1 kg/s, matching
the existing profile integrator. Zero points retains prescribed base-mass events.

For each event using a profile, supply `injection_end_times(e)`. The code
integrates the linear profile from `injection_times(e)` through that end time,
then multiplies by `injection_scales(e)`. Both endpoints must lie within the
profile, with end >= start. This replaces the port's base event mass; its base
mass remains the template normalization. Integrated mass is bounded to 1e6 kg
and at most 1e6 times that normalization. Drop size/composition stay fixed and
multiplicity changes. Zero integrated mass consumes an event but no particle IDs.

The entire interval's mass is injected at its start time, including an event
starting at final time. This is discrete packet injection; choose event spacing
for the desired temporal resolution. Intervals may overlap if explicitly
specified. A profile does not create events automatically. Profiles and interval
ends are part of the restart context. The arrays are supplied in the namelist;
external rate-curve file loading is not provided here.

## EB wall impacts and stationary films

The default `wall_impacts=.false.` retains specular reflection. To enable
mixture deposition/splash and film heating/evaporation, add to `/eb_spray/`:

```fortran
wall_impacts=.true.,
particle_options%liquid_viscosity=.001,
particle_options%surface_tension=.072,
walls%enabled=6*.true.,walls%temperature=6*400.,
walls%evolve=.true.,walls%liquid_conductivity=.6,
```

Viscosity and surface tension must be positive. All six wall-temperature entries
must match: they represent one uniform EB wall temperature. Face-enable flags
do not select individual EB facets. Existing contact angle, film capacity,
minimum height and per-material conductivity members of `spray_wall_options`
remain available. `walls%evolve` requires impacts. Contact depth is limited to
32 per particle source substep. Film thermal evolution occurs in separate half
steps; films created in the final particle half begin evolving on a later step.

`result.csv.films.rank<R>.csv` contains material/port, position, normal, total
mass, temperature, diameter, height and total component masses. Empty ranks
write headers. Films have no global ID or canonical ordering. Files are protected
from overwrite. The budget includes film mass, species and energy plus the wall
impulse/heat ledger. Films follow the current finest fluid-side owners through
regrid and load redistribution.

Impact checkpoints use the V2 particle/film payload and retain splash lineage,
retired-parent ID watermark and film composition. The normalized wall settings
are part of the restart context. Checkpoints predating these new namelist fields
are rejected by strict context comparison. Continued runs can change rank count.

The focused fixture covers deposition at 20 m/s and splash at 200 m/s, using
1e-16 kg total injected mass for the splash case to remain in the bounded thin-film
range. The pinned wet-impact correlation can produce inadmissible child fractions
outside its supported range; this terminates with an explicit model-range error
instead of repeated timestep retries. This integration does not qualify general
thick-film impact, lateral film flow, nozzle or diesel combustion physics.
See [dated evidence](../../docs/validation/spray-mpi-eb-impact-application-20260918.md).
