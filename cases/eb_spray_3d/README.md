# Experimental serial EB spray application

Build with `PELEF_ENABLE_SPRAY=ON`, then run in a copy of this case directory:

```sh
/path/to/build/pelef_eb_spray_3d gas.nml spray.nml
```

The fixed executable uses the elementary H2/O2 mechanism and its existing default
mole composition. For a configured mechanism, the selected EB executable's name
has `_spray` appended (normally `pelef_reactive_eb_3d_selected_spray`). Set
`thermo_model='selected'` and the existing named composition fields in `gas.nml`.
Both frontends share the application driver.

`gas.nml` uses the existing `reactive_eb_3d` schema. This application requires
`geometry_type='level_set'`: the builder returns the original nodal field to the
particle solver, binding gas and particle geometry. `level_set.dat` begins with
cell counts and then nodal values in Fortran order (x fastest); positive is fluid.
The supplied plane is x=1.25 mm in a 4 mm cube with 1 mm cells.

`spray.nml` has an `eb_spray` group:

- `nc` (1–16), `vapor_names`, `densities`, `heat_capacities`,
  `reference_temperatures`, `reference_pressures`, `latent_heats` and `fractions`
  define synthetic constant-property liquid components. Active fractions sum
  to one; inactive fractions are zero. Vapor names must exist in the mechanism.
- `number`, `total_mass`, `center`, `axis`, `speed`, `half_angle`, `diameter` and
  `liquid_temperature` configure each cone injection. `number` is
  limited to 125000. Mixture density sets per-droplet mass; multiplicity sets
  represented injected mass.
- `injection_events` (1–1024), `injection_times` and `injection_scales` define
  nondecreasing event times in seconds and nonnegative mass multipliers. Defaults
  retain one event at time zero with multiplier one. Ties are processed in array
  order. Zero multipliers advance the event counter without consuming IDs. Events
  at the final time are included; later events remain unprocessed. Each positive
  event injects `number` parcels representing `total_mass * injection_scales`.
  Every candidate position must be in fluid before population/IDs/inventory commit.
  The total retained population must fit both 125000 and `options%maximum_parcels`.
- `options` is the existing `spray_options` derived namelist object. The API
  supports representative-drop breakup; KH child creation is rejected.
- `maximum_dt` limits the gas step in addition to existing EB hydro and transport
  bounds. Failed coupled steps retry at most 30 halvings. Particle source/drift
  substeps follow their own limits.

The sample injects two water parcels 1 nm from the wall, toward it at 20 m/s.
The 20 ns run exercises reflection, heating and evaporation with gas chemistry
and molecular transport. It writes the existing gas CSV to `output_file` and
particle/component fields to `<output_file>.particles.csv`. Standard output
reports accepted steps/time, processed events and cumulative injected, exterior,
wall and chemistry inventories. CSVs are analysis outputs, not restart files. Output paths must
not alias the gas, spray or geometry inputs; checks normalize lexical absolute
paths, not filesystem symlinks.

## External rate files

Rate tables can also come from a separate optional group in `spray.nml`:

```fortran
&eb_spray_rate_files
  profile_file='shared.dat',profile_files(2)='port2.dat'
/
```

`profile_file` supplies the shared profile; `profile_files(p)` overrides it for
port `p`. Files contain 2--64 whitespace-separated pairs of time in seconds and
rate in kg/s, with blank lines and `#`/`!` comments allowed. Times must strictly
increase within 0--1e6 seconds; rates must be within 0--1 kg/s. Leave the
corresponding inline tables empty and counts at their defaults (port count zero
is also accepted). Conflicting inline/file data and unused port paths reject.
Existing event intervals and scaling still apply; files do not create events.

Paths are relative to the process working directory. Files load at startup,
including restart. Loaded numbers enter checkpoint identity; file paths and
comments do not. The shared file follows the same inheritance rules as the
shared inline profile. See [focused checks](../../docs/validation/spray-serial-eb-rate-files-20260925.md).

## Save and resume

The existing gas settings `checkpoint_interval_steps=N` and
`checkpoint_file='run.chk'` publish immutable snapshots named `run.chk.stepN`,
`run.chk.step2N`, etc. `stop_after_checkpoint=.true.` exits after the first new
snapshot. Saving happens after all injections at the accepted time, before final
CSV output. Resume with `restart_file='run.chk.stepN'`; the saved step is not
immediately saved again. Use a fresh prefix when rerunning a save, because an
existing snapshot is never overwritten.

Snapshots include gas, temperature, every parcel field, component masses, next
event, last ID, initial gas inventory and cumulative injected/exterior/wall/
chemistry ledgers. Restart skips initial injection and resumes pending events.
Final time, maximum step count, output paths and checkpoint controls may change.
Geometry values, gas configuration, mechanism/transport tables, liquid options,
spray schedule and timestep controls must match. The context currently has a
131072-character capacity; large geometry/mechanism combinations can exceed it
and are explicitly rejected. This is a bounded serial format, with Linux atomic
publication inherited from the spray checkpoint API, not a power-loss durability
or cross-compiler restart guarantee. Covered cells retain valid EOS states from
the initializer and are checked by the existing checkpoint reader.

### Optional EB wall impacts

Add these settings to `&eb_spray` to replace specular reflection with the
existing mixture impact model:

```fortran
 wall_impacts=.true.
 options%surface_tension=.072, options%liquid_viscosity=.001
 walls%temperature=300, walls%maximum_films=1000
 walls%evolve=.true., walls%conductivity=.6, walls%minimum_height=1e-7
```

These are illustrative water-like constants, not calibrated properties.
`walls%maximum_hits` defaults to 32 per trajectory. Films are written to
`<output_file>.films.csv` (header only when empty). Its component masses are
total kg; parcel component masses remain kg per represented drop. Snapshots
now also retain films and their components, plus retired-parent ancestry for
splashed children. Context V2 binds the impact settings and rejects older V1
snapshots. New deposits begin film exchange on the next source substep.

The inherited wet-impact correlations have a restricted valid range. For
example, normal incidence at film thickness/drop diameter 0.1123 produces a
negative child mass fraction and is rejected. The application now reports
`EB impact correlation outside admissible mass/velocity range` immediately,
without timestep retries for that model rejection. This remains an open model limitation, not
qualified general wet-wall impact support. See
[application evidence](../../docs/validation/spray-eb-impact-application-20260918.md).

Particle outflow is opt-in as described below. Other physical outer-boundary
particle modes and distributed EB impact/film ownership remain separate work.
Without particle outflow, trial trajectories must stay in the root domain. The
particle source/drift split is first order; this synthetic implementation case
is not calibrated spray validation or a convergence study.

### Computational parcel splitting

With ETAB or representative-drop KH-RT, set `options%maximum_multiplicity` to a
positive cap. The application splits parcels after each complete reactive update;
zero disables splitting. `options%maximum_parcels` bounds the resulting population
(also limited to 125000 here). Children copy per-drop component masses and retain
parent lineage; represented mass, momentum, energy and extensive KH state are
partitioned conservatively. IDs start above the shared watermark, including IDs
of retired impact parents. Capacity or ID overflow rejects the whole step.

Splitting-enabled checkpoints append `_SCOUNT_V1` to their format identifier and
store cumulative split IDs immediately after the ID watermark. Restore subtracts
that count before checking injection and four-child splash groups. Disabled-mode
checkpoint layout is unchanged. Older splitting-enabled snapshots without this
counter are rejected. Final-time injections are retained unsplit until the next
reactive update. See [focused evidence](../../docs/validation/spray-eb-split-20260918.md).

### KH shed-mass children

Set `options%breakup_model='khrt'` and `options%khrt_children=.true.` to
create KH shed-mass children. The shared watermark covers injection, KH,
splash and computational splitting. Children and pending KH inventories survive
checkpoint/restart. Restore checks four-child release groups after subtracting
injected and computational split IDs. The existing checkpoint layout already
stores the KH fields; no new format is needed for KH mode alone.

A crossing parent skips KH for that source substep and takes the wall path.
Other parents and their children drift after release; the transverse child energy
is balanced against the gas. This is source-first splitting, not a calibrated
breakup validation. See [focused evidence](../../docs/validation/spray-eb-khrt-application-20260918.md).

### Finite nozzle and size distributions

The shared scheduled nozzle accepts `nozzle_radius` in metres (default zero),
`nozzle%profile='solid'` or `'hollow'`, and `nozzle%swirl` and `nozzle%spread`.
These use the existing deterministic cone sampler. For example:

```fortran
 nozzle_radius=1e-4, half_angle=.1
 nozzle%profile='hollow', nozzle%swirl=.4, nozzle%spread=.02
 sizes%distribution='lognormal', sizes%std=2e-6
 sizes%minimum=5e-6, sizes%maximum=2e-5, sizes%weighting='number'
```

`sizes%shape` enables the existing Weibull sampler; the other existing choices
are `normal`, `lognormal` and `chi-squared`. The default remains monodisperse.
`sizes%weighting='mass'` assigns equal represented mass to records; `'number'`
assigns equal represented drop counts. The requested event dose is preserved.
Every generated point must lie in fluid; validation precedes population and
budget publication. All controls enter normalized restart context, so checkpoints
from before these added fields fail context comparison. This does not add
independent nozzle ports. See [focused evidence](../../docs/validation/spray-eb-nozzle-20260918.md).

### Component property models

The component table accepts the same property controls as the MPI EB application:

| Model | Inputs |
| --- | --- |
| Density | `mixture_density_models`, `mixture_density_coefficients(4,16)`, `mixture_density_temperature_ranges(2,16)` |
| Heat capacity | `mixture_cp_models`, `mixture_cp_coefficients(4,16)`, `mixture_cp_temperature_ranges(2,16)` |
| Vapor pressure | `mixture_vapor_pressure_models`, `mixture_antoine_coefficients(4,16)`, `mixture_antoine_temperature_ranges(2,16)` |
| Vapor surface cap | `mixture_boiling_models`, `mixture_critical_temperatures` |

Density and heat capacity default to `constant`; `polynomial` uses coefficients
in ascending powers of temperature in kelvin and requires a valid temperature
interval. Vapor pressure defaults to `clausius`; `antoine` uses A, B, C, D for
`p_sat = D * 10**(A - B/(T+C))`, with pressure in pascals. The existing
`surface-cap` boiling option requires a critical temperature and represents a
bounded vapor-state cap, not general mixture boiling.

Unused component entries must retain defaults. Injection uses the mixture's
density, heat capacity and energy at its specified temperature; evolving parcels
and films use the component models. All property inputs enter restart context.
Older checkpoints lacking these fields are rejected. The supplied coefficients
must describe the intended liquid; illustrative test fits are not calibrated
material data. See [focused evidence](../../docs/validation/spray-eb-properties-20260918.md).

### Injection-rate histories

Set `injection_profile_points` to 2–64 and supply increasing
`injection_profile_times` (s) and nonnegative `injection_profile_rates` (kg/s).
Each event integrates the piecewise-linear profile from `injection_times(i)` to
`injection_end_times(i)` and multiplies the result by `injection_scales(i)`.
Injection still occurs as a batch at the event start; the integration interval
specifies its dose, not continuous particle insertion. All intervals must lie
within the profile. Zero dose consumes no IDs. With zero profile points, the
existing `total_mass * injection_scales(i)` rule applies.

For a finite nozzle, `velocity_from_rate=.true.` computes speed as scaled event
mass divided by interval duration, mixture density, circular outlet area and
`discharge_coefficient`. The coefficient defaults to 1 and must be in [1e-6,1];
rate-driven radius must be at least 1e-12 m. Positive event speeds cannot exceed
1e7 m/s. Prescribed `speed` remains the default when rate-driven velocity is off.
Example:

```fortran
 injection_profile_points=2
 injection_profile_times=0.,4e-8, injection_profile_rates=.001,.003
 injection_times(1)=0., injection_end_times(1)=5e-9
 nozzle_radius=1e-4, velocity_from_rate=.true., discharge_coefficient=.8
```

Unused profile slots remain zero; without a profile, end times remain -1. Profile
and velocity controls enter restart context, so older snapshots lacking those
fields are rejected. See [focused evidence](../../docs/validation/spray-eb-rate-20260918.md).

### Multiple nozzle ports

Set `port_count` to 1–16 and choose `event_ports(i)` for each scheduled event.
Unspecified ports inherit the shared injection settings. Per-port overrides are:
`port_numbers`, `port_centers(3,16)`, `port_axes(3,16)`, `port_masses`,
`port_speeds`, `port_angles`, `port_diameters`, `port_temperatures` and
`port_fractions(16,16)`. Fractions use the shared component table and must sum to
one for each active port. Defaults use -1 for unspecified counts/scalars/fractions
and the internal large-value sentinel for position/direction vectors.

Nozzle shape, size distribution, radius, rate profile and discharge coefficient
inherit the shared values unless overridden below.
With a rate profile, its integrated dose replaces
`port_masses`; rate-driven velocity uses each selected port's mixed density.
Every parcel carries its port in `material_id`, now also included in particle CSV.
Children and films retain that identity. Checkpoint restore validates port IDs
and counts only positive-dose events using the selected port's population.
All resolved port settings enter restart context; earlier contexts are rejected.
See [focused evidence](../../docs/validation/spray-eb-ports-20260918.md).

### Per-port rate profiles and outlets

`port_profile_points(16)` defaults to -1, meaning inherit the common rate profile.
Set a port to zero to use its prescribed event mass, or 2–64 to use
`port_profile_times(64,16)` and `port_profile_rates(64,16)`. Explicit profile
columns require their own point count; unused entries remain zero. Every event
uses its selected port's profile and the existing `injection_end_times(i)`.
A prescribed-mass event retains end time -1.

`port_radii(16)` and `port_discharge_coefficients(16)` also default to -1 and
inherit their shared values. Radius feeds both parcel positions and rate-derived
speed. `velocity_from_rate` remains a common switch: when enabled, every active
port requires a profile, positive radius and valid discharge coefficient. All
resolved profiles and outlet settings enter restart context. Example overrides:

```fortran
 port_profile_points(2)=2
 port_profile_times(1:2,2)=0.,4e-8, port_profile_rates(1:2,2)=.001,.003
 port_radii(2)=2e-4, port_discharge_coefficients(2)=.5
```

See [focused evidence](../../docs/validation/spray-eb-port-profiles-20260918.md).

### Per-port nozzle and size models

`port_nozzles(16)` and `port_sizes(16)` accept the same fields as `nozzle` and
`sizes`. Each omitted field inherits its shared counterpart. For example, a
second hollow/lognormal port beside a shared solid/Weibull port can specify:

```fortran
 port_nozzles(2)%profile='hollow', port_nozzles(2)%spread=.02
 port_sizes(2)%distribution='lognormal', port_sizes(2)%weighting='number'
 port_sizes(2)%std=2e-6, port_sizes(2)%shape=0.
```

Set inherited distribution-specific parameters explicitly when changing models
(e.g. `std=0` for Weibull, or `shape=0` for lognormal). The internal inheritance
markers are `'inherit'` for strings, -1 for nonnegative numeric controls and
`huge(real)` for signed swirl. Resolved values enter restart context; snapshots
from before these fields fail context comparison. See
[focused evidence](../../docs/validation/spray-eb-port-models-20260918.md).

### Particle outflow

`particle_outflow=.true.` allows parcels to leave through all six exterior faces,
matching the existing extrapolated-outflow gas boundary convention. EB walls
use specular reflection by default, or the impact/film model when
`wall_impacts=.true.`. The default is false and retains root-contained trajectories.

Without impact films, escaped records retain their IDs with zero mass and zero
component inventory. With impacts enabled, escaped records are removed alongside
consumed parents; the shared ID watermark remains unchanged by removal.
Mass, momentum, liquid internal plus kinetic energy and mapped species losses
join `Boundary inventory`, separately from `Wall inventory`. KH children and
computational splitting remain supported. The reactive update and its retries
publish all state and ledgers together. Existing checkpoint fields persist the
combined boundary history and dead records; the new mode flag enters restart
context. Earlier contexts are rejected. See
[specular evidence](../../docs/validation/spray-eb-outflow-application-20260918.md)
and [impact/film restart evidence](../../docs/validation/spray-eb-splash-outflow-application-20260918.md).

### External acceleration

Set `acceleration=0,0,-9.81` in `eb_spray` for a uniform acceleration in m/s².
The default is zero; each component must be finite and at most 1e6 in magnitude.
Gas and live particles receive impulse/work; stationary attached films do not.
`Force inventory` reports cumulative external exchange separately from wall,
exterior and chemistry budgets. Nonzero acceleration adds an `_EXT_V1`
checkpoint inventory and the setting is part of restart context. See
[restart evidence](../../docs/validation/spray-eb-acceleration-application-20260918.md).

## Artificial viscosity

The gas namelist accepts finite nonnegative `difmag` (default 0.1). Positive
values apply the EB gas flux correction between particle half steps and must
match when restarting. Negative values, including -1, are invalid. See the
[focused evidence](../../docs/validation/eb-spray-difmag-20260925.md) for checked
configurations and checkpoint-format compatibility limits.

The EB coefficient default changed to 0.1; old zero-policy checkpoints require
explicit `difmag=0`. See [default-policy checks](../../docs/validation/eb-default-difmag-20260925.md).
