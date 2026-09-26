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

## Recursive SDC

The MPI spray executable selects recursive gas/particle stepping with these
`/eb_tree/` controls:

```fortran
time_integrator='sdc', sdc_iterations=2, subcycling=.true.,
hydro_method='plm_ctu', riemann_solver='pelec', boundary_faces=6*'outflow',
```

This experimental path supports drag, heating and evaporation with gas chemistry
and molecular transport. Gas force density uses `external_forcing`;
`acceleration` in `/eb_spray/` accelerates both gas and particles. Particle
half kicks bracket each level interval, with physical-particle momentum and work
reported in the `force` budget column
([checks](../../docs/validation/spray-recursive-acceleration-20260927.md)). Keep
`particle_options%breakup_model='none'`, `wall_impacts=.false.`,
`walls%evolve=.false.` in `/eb_spray/`; recursive breakup and films reject at input. Periodic migration also rejects.
Reflecting particle/EB walls are distinct from the unsupported film-impact model.

Physical particle exchange is synchronized conservatively after each outer step.
This is a PeleF extension to the pinned source-transfer algorithm. See the
[application checks](../../docs/validation/spray-recursive-app-20260927.md) for the
bounded test scope. Recursive checkpoints carry a separate integration-mode tag;
a focused two-level case now checks regridding followed by 2-to-4-rank
continuation with a pending injection event
([evidence](../../docs/validation/spray-recursive-restart-regrid-20260927.md)).
Three-level evolution and broader restart/model combinations remain open.

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
regrid controls. Optional KH shedding and wall impacts create child parcels. It inherits root-contained EB flight and the
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
ordered table; each port can also load a separate rate-curve file.

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
ends are part of the restart context.

Alternatively, append a separate optional group to `spray.nml`:

```fortran
&eb_spray_rate_files
  profile_files(1)='rate1.dat',profile_files(2)='rate2.dat'
/
```

Each file contains 2--64 whitespace-separated `time_seconds rate_kg_per_second`
pairs. Blank lines and `#`/`!` comments are allowed. The same profile bounds
apply. Leave that port's inline point count and tables at zero; supplying both
forms is rejected. Paths for unused ports are rejected. Paths are relative to
the process working directory; all MPI ranks must see consistent files.
Files load at startup, including restart. Loaded numbers enter the existing
checkpoint context, so changed values reject continuation while changed paths
or comments do not. Event times and interval ends must still be supplied.
See [focused checks](../../docs/validation/spray-eb-rate-files-20260925.md).

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

## Finite outlets and size distributions

Each port accepts `port_radii(p)` in metres, `port_nozzles(p)` and
`port_sizes(p)`. Defaults retain a point outlet and one prescribed diameter.
For example:

```fortran
port_radii(1)=1.e-4,
port_nozzles(1)%profile='solid',port_nozzles(1)%swirl=.4,
port_sizes(1)%distribution='weibull',port_sizes(1)%shape=2.5,
port_sizes(1)%minimum=5.e-6,port_sizes(1)%maximum=2.e-5,
```

Nozzle profiles are `legacy`, `solid` and `hollow`. Solid/hollow profiles require
a positive radius. `swirl` is a rotation in radians of the transverse velocity
relative to the outlet radial direction. Hollow outlets place parcels at the
outlet radius; optional `spread` sets their angular spread about the port cone
angle. Solid outlets use the existing deterministic radial/cone quadrature.
The supplied finite outlet must lie inside fluid; EB migration rejects samples
outside the domain or inside the solid.

Size fields are `distribution` (`weibull`, `normal`, `lognormal`, `chi-squared`),
`weighting` (`mass` or `number`), `shape`, `std`, `minimum` and `maximum`.
Use `shape` for Weibull and `std` for normal/lognormal; unused parameters must
satisfy the existing injection-kernel contract. Diameters and bounds use metres;
`port_diameters(p)` remains the distribution's nominal diameter. The default
Weibull shape zero retains the single prescribed size. Mass weighting gives
records equal represented mass; number weighting gives equal multiplicity.
Both normalize to the port's event mass before profile/scaling factors apply.

The port template is deterministic and reused at successive events, with new
global IDs. These are injection geometry and size controls, not stochastic
particle dispersion. They are included in the restart context; checkpoints
predating the new fields fail the strict context comparison. See
[focused evidence](../../docs/validation/spray-mpi-eb-nozzle-20260918.md).

## Speed from the port mass-flow profile

Set `port_velocity_from_rate(p)=.true.` for a port with a positive finite radius
and a mass-flow table. `port_discharge_coefficients(p)` defaults to one and must
be in [1e-6, 1]. Other ports retain their prescribed speeds. A nondefault
coefficient on a prescribed-speed port is rejected.

For each event, the code integrates that port's rate table, applies the event
scale, and computes

```text
speed = event_mass / (event_duration * liquid_density * pi * radius**2 * discharge_coefficient)
```

Density is evaluated from the port's liquid composition and temperature.
`event_duration` is `injection_end_times(e)-injection_times(e)`. The speed uses
the average scaled rate over this interval, while the whole packet is injected
at the interval start, as for the existing profile-mass input. The port's
prescribed speed is overridden; cone and swirl directions are retained.
The resulting velocity enters both particles and injected momentum/energy.

A zero-mass event reserves no IDs and injects nothing. Positive events require
a positive duration and speed at most 1e7 m/s, the existing injection bound.
Radius and coefficients are validated before stepping. The mode, coefficients,
profiles and event intervals are bound into restart context. Example:

```fortran
port_velocity_from_rate(1)=.true.,port_discharge_coefficients(1)=.8,
port_radii(1)=1.e-4,
port_profile_points(1)=2,
port_profile_times(1:2,1)=0.,4.e-8,port_profile_rates(1:2,1)=.002,.002,
injection_times(1)=0.,injection_end_times(1)=5.e-9,
```

See [focused evidence](../../docs/validation/spray-mpi-eb-rate-velocity-20260918.md).

## Component property models

The shared `nc` component table accepts the same `mixture_*` model controls as
the regular spray application:

| Model | Fields |
| --- | --- |
| Density | `mixture_density_models`, `mixture_density_coefficients(4,16)`, `mixture_density_temperature_ranges(2,16)` |
| Heat capacity | `mixture_cp_models`, `mixture_cp_coefficients(4,16)`, `mixture_cp_temperature_ranges(2,16)` |
| Saturation pressure | `mixture_vapor_pressure_models`, `mixture_antoine_coefficients(4,16)`, `mixture_antoine_temperature_ranges(2,16)` |
| Vapor surface cap | `mixture_boiling_models`, `mixture_critical_temperatures` |

Density and heat capacity default to `constant`; `polynomial` selects four
coefficients in ascending powers of absolute temperature. Supply a valid fit
interval. Existing validators reject nonpositive fits and out-of-range states.
The base component reference enthalpy still comes from its gas species and the
supplied latent heat; the polynomial heat capacity determines energy changes.

Vapor pressure defaults to `clausius`. `antoine` uses coefficients A, B, C, D
for `p_sat = D * 10**(A - B/(T+C))`, with T in kelvin and D in pascals.
`mixture_boiling_models='surface-cap'` uses the existing Clausius surface-state
cap and requires a valid critical temperature; it is not a general physical
mixture-boiling model. Unused component model entries must retain defaults.

Each port's composition and temperature determine its mixed density, heat
capacity and energy. These values set injected drop mass and rate-derived speed;
the component models also feed evaporation, heating, breakup, impacts and film
exchange through the existing mixture APIs. All model parameters enter restart
context. Checkpoints preceding these fields fail strict context comparison.
See [focused evidence](../../docs/validation/spray-mpi-eb-properties-20260918.md).

## KH-RT child parcels

The MPI EB executable now accepts the existing KH child controls, for example:

```fortran
particle_options%breakup_model='khrt',particle_options%khrt_children=.true.,
particle_options%liquid_viscosity=.001,particle_options%surface_tension=.072,
```

The shared ID watermark covers both injection and four-child releases. Child
composition, lineage and pending KH state follow migration, regrid and checkpoint
restore. The transverse kinetic-energy contribution is balanced against the gas.
A source parent crossing the EB skips KH release for that substep and takes the
wall path; this applies to both specular reflection and optional impact films.
Unhit parents and generated children advance through the remaining source flight.
The source operator uses breakup followed by drift, not globally ordered physical
breakup times. The full gas/particle/film update retains rollback, including IDs.

Checkpoint restore verifies the watermark against processed injections plus
four-child release groups. It works with the V1 particle-only payload and V2
impact/film payload. The mode remains experimental and its constants need
physical calibration. See [focused evidence](../../docs/validation/spray-mpi-eb-khrt-application-20260918.md).

## Computational parcel splitting

For ETAB or KH-RT, `particle_options%maximum_multiplicity>0` enables splitting
after both particle source halves of each accepted reactive step. Zero disables
it. Descendants retain their final owner, composition and lineage; represented
mass, momentum and energy are preserved. `particle_options%maximum_parcels`
bounds each rank's population. Capacity failure rolls back the entire step,
including gas, films, particle IDs and budgets.

The checkpoint gas prefix additionally stores the cumulative number of IDs
created by computational splitting. Restore subtracts this count before checking
the injection and four-child KH/splash ID groups. The particle payload remains
V1 without impact films and V2 with them. Splitting-disabled checkpoints retain
their two-field injection metadata. Newly injected parcels are split by the next
reactive update; injection at the final time is not immediately split.
See [focused evidence](../../docs/validation/spray-mpi-eb-split-20260918.md).

## Particle outflow

Set `particle_outflow=.true.` in `eb_spray` to let particles leave through all
root faces. The default is false. This works with default specular EB reflection
and with `wall_impacts=.true.`. Removed particles do not recycle IDs; gas boundary
budgets include their mass, momentum, energy and vapor-species inventories while
EB wall budgets remain separate. The mode enters restart context. An empty
particle population can be checkpointed and resumed with a different rank count.
See [specular restart evidence](../../docs/validation/spray-mpi-eb-specular-application-20260918.md)
and [impact-film restart evidence](../../docs/validation/spray-mpi-eb-outflow-application-20260918.md).

## External acceleration

`acceleration=0,0,-9.81` in `eb_spray` applies uniform acceleration in m/s² to
gas and live particles. Each component must be finite and bounded by 1e6;
the default is zero. Attached films remain stationary. The budget CSV has a
`force` column, separate from wall and exterior exchange. Nonzero acceleration
adds force history to the coupled gas checkpoint prefix; its serialized input
binds the mode to restart. See
[rank-changing restart evidence](../../docs/validation/spray-mpi-eb-acceleration-application-20260918.md).

## Artificial viscosity

`&eb_tree` accepts finite nonnegative `difmag` (default 0.1). The gas update
uses it between particle half steps, and restart requires the same value.
The library requires coefficient agreement across ranks before particle work.
See [focused checks](../../docs/validation/mpi-eb-spray-difmag-20260925.md).

The EB coefficient default changed to 0.1; old zero-policy checkpoints require
explicit `difmag=0`. See [default-policy checks](../../docs/validation/eb-default-difmag-20260925.md).

### CTU gas update

The tree input accepts `hydro_method='plm_ctu'` with `riemann_solver='pelec'`.
The coupled MPI gas update supports common or recursive intervals. Method
identity is checked when restoring the gas snapshot. Focused fixed-mechanism
2-to-4-rank continuation and coupled conservation checks are recorded in
[the CTU evidence](../../docs/validation/mpi-eb-spray-ctu-20260925.md); these do
not establish nozzle, diesel ignition or validated LES behavior.

## Prescribed gas inflow

The tree gas input accepts explicit `boundary_faces`, `inflow_density`,
`inflow_pressure` and `inflow_velocity`, as described in the
[tree case guide](../reactive_eb_tree_3d/README.md#prescribed-inflow-and-outflow).
Set `particle_outflow=.true.` in `spray.nml`. Both gas inflow and outflow faces
allow existing particles to leave; the gas inlet does not create particles.
Changed gas boundary controls reject restart. See the
[fixed-mechanism MPI checks](../../docs/validation/mpi-eb-spray-inflow-20260926.md).

## Periodic boundaries

Set paired `periodic` labels in the tree input and `particle_outflow=.true.`
in the spray input. For the supplied y-normal plane geometry:

```fortran
 boundary_faces='periodic','periodic','outflow','outflow','periodic','periodic'
```

Gas and particle exchange use wrapped donors. Parcels crossing periodic faces
wrap into the domain and migrate to the resulting gas-patch owner. Other
explicit faces allow particle outflow. Restart retains the boundary labels;
changing them is rejected. `wall_impacts=.true.` also works with periodic faces:
flight resumes on the opposite face, EB impacts query periodic film images, and
film heating/evaporation samples and deposits to wrapped gas donors. Periodic
crossings add no exterior impulse or lost mass; physical outflow still does.
Opposing nodal geometry values must match. Stationary film positions remain
inside the root domain and migrate with regridding; their footprint can cross
the periodic seam.
KH child generation is supported. Parent eligibility checks the complete wrapped
flight with the same drift routine used for movement. Periodic crossings permit
release; EB or reflecting domain-wall contacts suppress release for that substep.
Children then drift and migrate with their parent population. Specular reflection
from the embedded surface remains available.

## Adiabatic domain walls

The tree input accepts `slip_wall` and `no_slip_wall` gas labels, for example:

```fortran
 boundary_faces='no_slip_wall','no_slip_wall','outflow','outflow','no_slip_wall','no_slip_wall'
```

Set `particle_outflow=.true.` in the spray input to enable explicit particle
boundary routing. Particles reflect specularly at either wall label; the gas
uses the corresponding slip/no-slip condition. Exterior particle impulse is
included in the coupled boundary inventory. Wall labels are bound to restart.
`wall_impacts=.true.` enables deposition/splash and stationary films on the
embedded surface while domain faces remain specular. Outer-face impulse enters
the boundary budget, and EB impact/film exchange enters the wall budget. Flight
continues through both event types in time order. This does not create films
on the outer domain faces. KH child generation
is supported for nonperiodic faces: a parent reaching a domain wall skips KH
release on that substep, reflects, and can release on later eligible substeps.
Child state and ID allocation participate in the existing interval transaction.

`isothermal_wall` is also accepted with `wall_temperature(6)` in the tree input:

```fortran
 boundary_faces='isothermal_wall','isothermal_wall','outflow','outflow','isothermal_wall','isothermal_wall'
 wall_temperature=6*2000
 transport_enabled=.true.
 automatic_timestep=.true.
```

Supply six finite positive temperatures in kelvin in x-low/high, y-low/high,
z-low/high order. Gas molecular conduction uses the isothermal face values;
particles reflect specularly without a direct particle-wall heat model.
Gas heat transfer and particle domain-wall impulse enter the boundary ledger.
The common/recursive coupled transaction includes both processes. Restart binds
all six wall temperatures; changed values reject before continuation. The
domain-wall film restriction above also applies to thermal walls.

## Common-step SDC

Set `time_integrator='sdc',sdc_iterations=2,hydro_method='plm_ctu',
riemann_solver='pelec'` in `&eb_tree`, with explicit `boundary_faces`.
Set `particle_outflow=.true.` in `&eb_spray`. This connects gas SDC,
particle half kicks/full drift, internal source exchange and final computational
splitting. `external_forcing` is gas-only force per volume (N/m³); spray
`acceleration` is acceleration (m/s²) applied to both phases. The gas momentum
source is `external_forcing + density*acceleration`; its work uses the local
gas velocity.

Coupled checkpoints preserve reaction history and total external impulse/work
alongside parcel state, IDs and injection progress. Common-step SDC currently
requires `subcycling=.false.`. `wall_impacts=.true.` now connects EB
deposition/splash and film evolution using the existing `walls` controls; film
state and splash IDs are retained across restart. Recursive film/splash scheduling
remains unimplemented; recursive drag/heating/evaporation uses the controls above. The normal split integrator remains
the default. This experimental path is not physical spray validation.
