# Experimental MPI AMR spray application

Enable `PELEF_ENABLE_MPI=ON` and `PELEF_ENABLE_SPRAY=ON`. Run from a fresh output
folder, passing absolute input paths:

```sh
mpiexec -n 2 /path/to/build/pelef_mpi_amr_spray_3d \
  /path/to/cases/mpi_amr_spray_3d/gas.nml \
  /path/to/cases/mpi_amr_spray_3d/tree.dat \
  /path/to/cases/mpi_amr_spray_3d/spray.nml
```

The standard `reactive_3d` gas namelist controls the gas initial condition,
domain and per-face gas boundaries, chemistry/transport, numerical scheme, CFL and final time.
The topology file uses the existing tree format: child count followed by
`parent ilo ihi jlo jhi klo khi ratio` records, with parents before children.
Patch owners are assigned round-robin; ranks without gas patches are supported.

The separate `amr_spray` namelist supplies up to 16 liquid components by vapor
species name, mass fractions, density, heat capacity, reference temperature /
pressure and latent heat. `parcel_count` and `total_mass` describe one initial
cone injection on rank zero. `center`, `axis`, `speed`, `half_angle` (radians),
`diameter` and `drop_temperature` define the injected population. Its collective
ID watermark and component inventories follow particles through migration and
breakup. `maximum_parcels` bounds per-rank storage. Existing `spray_options`
fields can be set as namelist components, e.g. `options%breakup_model='etab'`,
with the required liquid viscosity/surface tension. `subcycling=.true.` selects
recursive gas steps and, when enabled, recursive LES half steps. The driver
retains a conservative finest-grid hydro/molecular-transport bound and selects
the scaled LES diffusion bound.

For timed injection, set `injection_events` (1 through 1024), nondecreasing
nonnegative `injection_times`, and matching `event_parcels` and `event_mass`
arrays. `event_ports` selects the port for each event (default 1). This
explicit schedule replaces the default time-zero injection; when
`injection_events=0` (default), `parcel_count` and `total_mass` retain their
original time-zero meaning. A zero-mass event must also have zero parcels.
For example:

```fortran
 injection_events=2, injection_times=0.,1.e-8
 event_parcels=4,4, event_mass=4.e-11,6.e-11
```

The solver lands exactly on each event. If regridding and injection coincide,
the new grid is established first. Equal-time events execute in array order.
Events at final time are applied before
output; later events remain unapplied. New parcels, component masses, IDs and
the cumulative injection ledger publish only after injection and migration
succeed. The summary reports the number of applied events.

Set `port_count` (1 through 16) to configure multiple outlets. `port_centers`
 and `port_axes` have shape `(3,16)`; `port_speeds`, `port_half_angles`,
`port_diameters`, `port_temperatures` and `port_radii` have one entry per port.
`port_fractions(1:component_count,port)` selects its composition from the shared
liquid table; omitted columns inherit `fractions`. For example:

```fortran
 port_count=2, event_ports=1,2
 port_centers(:,1)=.00075,.00075,.00075
 port_centers(:,2)=.00325,.00325,.00325
 port_axes(:,1)=1.,0.,0., port_axes(:,2)=-1.,0.,0.
 port_speeds=20.,30., port_diameters=1.e-5,2.e-5
 port_temperatures=350.,340.
 port_fractions(1:2,1)=.25,.75, port_fractions(1:2,2)=.75,.25
```

The existing nozzle model is selectable per port using
`port_nozzles(port)%profile` (`legacy`, `solid`, `hollow`), `%swirl` and `%spread`.
Positions on a finite outlet use `port_radii`. Explicit ports default to origin,
x direction, 10 m/s, zero cone angle/radius, 10 micrometer diameter and 350 K;
set their geometry explicitly. With `port_count=0` (default), the original
scalar cone controls populate port 1. Each newborn parcel stores its port in
`material_id`, which survives migration and restart; mixture properties continue
to use the component inventory and shared component tables. Port settings bind
the checkpoint context.

An optional single scheduled regrid is specified by `regrid_time` and
`regrid_topology` in `amr_spray`. The latter path is relative to the run directory
unless absolute. The solver lands on that time and atomically remaps gas and
particles. The root physical domain/resolution remains fixed. `root.dat` is a
coarsening target; `tree.dat` refines one corner.

For temperature-driven adaptation, set `adaptation_interval` to a positive
number of completed gas steps. Configure `refinement_levels` (0 through 16),
`refinement_ratios`, `temperature_thresholds`, `tile_shape`, `tag_buffer` and
`maximum_nodes`. Thresholds must be positive and nondecreasing; ratios are at
least two. For example:

```fortran
 adaptation_interval=1, refinement_levels=1, refinement_ratios=2
 temperature_thresholds=1005., tile_shape=1,1,1, tag_buffer=0, maximum_nodes=100
```

This uses the existing composite-temperature tagger. When topology changes,
gas and particles remap together and owners are assigned round-robin. The
initial tree comes from the topology input; checks begin after completed steps,
including the final step. A zero interval disables adaptation. Automatic
adaptation and `regrid_time` are mutually exclusive. Adaptation occurs before
injection and saving at a shared step/time. On restart, the saved topology is
restored and the saved step's adaptation check is not repeated. The next due
completed step resumes checking; adaptation controls bind the checkpoint context.
This does not perform particle-density tagging.

Set `rebalance_interval` to a positive number of accepted gas steps to sample
existing per-patch gas work timers on those steps. `rebalance_minimum_gain`
(default 0.1) is the required fractional reduction in estimated peak rank work.
The driver resets timing on each retry and uses only the accepted attempt.
Ownership changes commit gas and particle migration together; receive-capacity
failure preserves their prior ownership/state. Rebalancing runs after the step,
before the next adaptation/injection/checkpoint processing. It preserves the
geometry and does not interpolate gas fields.

The measured cost covers instrumented gas operators, not particle exchange,
particle models or all MPI communication. Timing-dependent owner maps are not
promised to repeat exactly between runs. Checkpoints remain rank-neutral;
restart starts with round-robin owners and resumes measurements on the next
due step. These controls bind the checkpoint context. A zero interval disables
measurement/rebalancing. This functionality is not scaling qualification.

Each step calls particle half / reaction-transport-hydro / particle half,
retrying rejected whole steps at half the interval. Gas CSV contains exposed
composite cells only. `<output_file>.particles.<rank>.csv` contains each rank's
particle IDs, lineage, material, position, velocity, mass, temperature,
multiplicity and component masses. These CSVs are analysis output, not restart
files. Standard output records step/time, last ID, regrid status and injected
mass/momentum/energy/species inventory.

For a configured mechanism bundle, `pelef_mpi_amr_spray_3d_selected` uses the
same driver. Set `thermo_model='selected'` and the standard named gas composition;
liquid vapor species must exist in that bundle. Bundle SHA and chemistry
integrator are printed by the selected frontend.

The fixture uses synthetic water/peroxide properties to exercise implementation.
It is not physical spray validation. This executable supports periodic, reservoir inflow/outflow and wall
boundaries, initial or scheduled multiport injection and one optional prescribed
regrid or temperature-driven adaptation, and gas-work-based rebalancing. Wall
films are opt-in through `walls` below; EB handling is not connected here. Gas-only AMR checkpoint/adaptation
namelists are not controls for this executable. Use the `amr_spray` controls above.

## Checkpoint and changed-rank continuation

Add a separate group to the spray input:

```fortran
&amr_spray_checkpoint
 checkpoint_file='saved_state', checkpoint_time=1.e-8
 stop_after_checkpoint=.true.
/
```

The driver lands on the requested time and applies coincident regrid/injection
operations before saving. Without `stop_after_checkpoint`, it continues normally.
A prefix is immutable: choose a new one for every save. The complete checkpoint
consists of `saved_state`, `saved_state.gas` and `saved_state.particles`; keep
these together with the same prefix. `saved_state.pending` is the retained
reservation/manifest staging file. A failed save may leave incomplete files,
but it does not publish the `saved_state` completion marker. Use a fresh prefix
after such a failure. Publication uses POSIX hard links on a shared filesystem;
this is atomic visibility, not a power-loss durability guarantee.

To resume, replace the checkpoint group with:

```fortran
&amr_spray_checkpoint
 restart_file='/path/to/saved_state'
/
```

The MPI rank count may change. Gas nodes are reassigned round-robin and saved
particles migrate to their finest owners. Saved state includes topology, gas
U/T, every parcel/component and breakup/lineage field, time/step count, ID
watermark, next injection event, regrid status and cumulative injected inventory.
Consumed events are not repeated; future events continue from their original
schedule. The prefix/context and paired generation must match and malformed
payloads reject before output.

Gas/transport/liquid tables, numerical options, injection schedule, topology
inputs and spray settings are bound into the checkpoint context. Continuation
limits and gas output path may change; checkpoint controls are excluded from the
context. `maximum_parcels` remains per-rank capacity, and checkpoint gathering
currently limits the global population to one million parcels on rank zero.
Gas payload I/O remains distributed. CSV analysis files are not restart files.

### Physical boundaries and budgets

Set `boundary_faces` in the gas namelist in x-low/x-high/y-low/y-high/z-low/z-high
order; unset entries inherit `boundary_condition`. For example:

```fortran
 boundary_condition='periodic'
 boundary_faces='slip_wall','outflow','periodic','periodic','periodic','periodic'
```

Fixed and selected gas inputs accept `inflow_mach(6)`: negative values retain
reservoir inflow; 0 through 100 select Mach inflow on an `inflow` face. Set a
positive `inflow_temperature` and zero `inflow_velocity` for active faces.
`inflow_pressure=0` uses interior pressure; a positive value prescribes Pa.
Composition comes from the current donor gas. Hydro/molecular timestep limits
use current distributed ghosts and the finest-level bounds; LES uses its own
Mach-aware selector. The checkpoint context binds the inlet controls.

Particles reflect specularly at gas wall types and leave at inflow/outflow
faces. Gas reservoir state uses `inflow_temperature`, `inflow_pressure`,
`inflow_velocity` and optional named inflow composition. Without that composition,
it uses the initial gas composition. Gas `wall_temperature` and `wall_velocity`
retain their existing thermal/tangential wall meanings. These controls do not
add droplet thermal-wall exchange or EB impact physics. Film formation and
evolution require the explicit `walls` controls below.

`Boundary inventory:` reports cumulative gas-plus-particle mass, three momenta,
total energy and species exchanges in conserved-variable order. Positive values
add inventory to the domain. It is separate from `Injected inventory:` and is
saved/restored with the gas checkpoint. Faces and reservoir/wall data bind the
checkpoint context. Place injection outlets inside the physical domain;
injection migration rejects exterior points on nonperiodic axes.

### Opt-in wall impact and films

Add `walls` settings to `&amr_spray`, and select a compatible wall face in the
gas namelist. For a low-x wall, a synthetic setup is:

```fortran
 walls%enabled(1)=.true.,walls%temperature=300.
 walls%evolve=.true.,walls%liquid_conductivity=.6,walls%minimum_height=1.e-10
 walls%maximum_films=100000
 options%surface_tension=.072,options%liquid_viscosity=.001
```

All wall faces are disabled by default. Enabled faces use the existing planar
impact model to deposit films or produce reflected/splashed parcels. Positive
liquid viscosity at each port temperature and a valid surface tension are
required. `walls%evolve` enables film heating and evaporation; it requires at
least one enabled wall. Gas wall temperatures and film wall temperatures are
separate existing model controls, so configure them consistently for the intended
case. These example properties are implementation fixtures, not calibrated data.

Films start empty and form through impacts. Coupled steps, prescribed/adaptive
regridding and gas-work rebalancing carry film state transactionally. Film settings
bind the checkpoint context; continuation restores film state and cumulative wall
exchange. Newly serialized settings can reject checkpoints from older application
builds even when defaults were used; cross-build compatibility is not promised.

Each rank writes `<output_file>.films.<rank>.csv` with material, position, normal,
mass, temperature, diameter, height and component masses, including a header when
that rank has no films. Film patches have no global ID and may change output rank
or row order after restart. These CSVs are analysis output; restart uses the
checkpoint's V2 particle/film payload. No lateral film transport or EB impact
geometry is implemented in this path.

### Pressure outlets

Set `outlet_pressure(face)` in the gas namelist on an `outflow` face to enable
the existing acoustic pressure ghost condition, for example
`outlet_pressure(2)=100000.` for the high-x face. Zero retains extrapolation;
positive pressures must not exceed `1e15 Pa`. Fixed and selected applications
pass this condition through the gas split, LES steps and LES interval estimate.
Pressure controls already belong to the gas checkpoint identity, so changed
values reject restart. Particles retain their outflow-removal behavior. A large
pressure jump can produce inadmissible ghost states and be rejected; this is not
a general nonreflecting boundary. See
[focused evidence](../../docs/validation/spray-amr-pressure-20260918.md).

### Static LES

Both fixed and selected executables accept an optional group in the spray input:

```fortran
&amr_spray_les
 les_model='vreman',les_cs=.17,les_ci=.1
 les_prt=.7,les_sct=.7,les_cfl=.1
/
```

Models are `smagorinsky` (default), `wale` and `vreman`. `les_cs` is the
model's eddy coefficient in `[0,2]`; `les_ci` controls isotropic stress in
`[0,1]`. Both default to zero, disabling LES. PrT/ScT default to 0.7 and must
lie in `[1e-6,1e6]`. The diffusion CFL defaults to 0.1 and must be in `(0,0.2]`.
These example coefficients are implementation fixtures, not calibrated values.

LES runs as two half steps around the gas split, using common stepping by
default or ratio-squared level intervals with `subcycling=.true.`. Its initial
diffusion estimate follows the same mode and limits the application interval;
subsequent particle feedback and nonlinear stages can still require the existing
retry/halving path. The subcycling flag and nondefault LES settings bind the checkpoint
context, so changing them on restart is rejected. Default settings retain the
previous context serialization. EB grids are not supported here.
See [application evidence](../../docs/validation/les-amr-application-20260918.md).

### Dynamic LES

Periodic, prescribed reservoir inflow/extrapolated outflow and slip/adiabatic/isothermal wall cases can select the dynamic path:

```fortran
&amr_spray_les les_model='dynamic',les_cfl=.1 /
&amr_spray_dynamic
 les_test_model='box',les_test_ratio=2,les_grid_ratio=1,les_harmonic_mean=.false.
/
```

The second group is optional with the shown defaults. Dynamic mode is active
even with `les_cs=les_ci=0`; its computed coefficients replace static
`les_cs/les_prt/les_sct/les_ci`. Filter models/ratios follow the regular Dynamic
LES kernel. Pressure-target outlets apply acoustic extrapolation through the
full donor halo; an inadmissible outer donor rejects the step. Nondefault dynamic
filter controls on a static model reject. Dynamic settings enter checkpoint identity, so changing
them rejects restart. Static/default context serialization is preserved.
The same particle/gas transaction, two LES half steps and retry logic apply.
The diffusion estimate includes padded donor cells and can be conservative.
Fixed/selected builds pass. Focused fixed-mechanism MPI parity and 2-to-4-rank
restart checks pass ([runtime checks](../../docs/validation/mpi-dynamic-amr-runtime-20260925.md)).
Dynamic LES with adaptive regridding/rebalancing and selected H2O2 continuation
also pass focused checks ([combined evidence](../../docs/validation/les-dynamic-adaptive-20260925.md)).
Detailed-fuel and physical validation remain outside this evidence.
Moving walls require `transport_enabled=.true.,viscosity_enabled=.true.` in the gas input.
See [wall integration checks](../../docs/validation/les-dynamic-amr-walls-20260925.md).
See [pressure-outlet checks](../../docs/validation/les-dynamic-amr-pressure-20260925.md).
See also [inflow/outflow checks](../../docs/validation/les-dynamic-open-20260925.md).

### Uniform acceleration

Add an optional group to the spray input to apply a uniform acceleration to
gas and free parcels (SI units, each component bounded by `1e6`):

```fortran
&amr_spray_force
 acceleration=0.,0.,-9.81
/
```

The default is zero. Attached wall films stay stationary. The application
prints cumulative external momentum and energy exchange as `Source inventory:`
in conserved-variable order, separately from injection and boundary inventories.
Nonzero acceleration binds restart identity and adds a `SPRAY_SOURCE_V1` record
to the particle checkpoint payload, with a matching gas-context marker.
The source inventory restores together with gas and parcels on a changed rank
count. Zero acceleration retains the previous context and payload format.
See [focused evidence](../../docs/validation/spray-amr-force-application-20260918.md).

### Component property fits

The optional `&amr_spray_properties` group configures each liquid component in
`vapor_names` order, independently of port composition. For example, a synthetic
first-component density and heat-capacity fit is:

```fortran
&amr_spray_properties
 mixture_density_models(1)='polynomial'
 mixture_density_coefficients(:,1)=1100.,-.2,0.,0.
 mixture_density_temperature_ranges(:,1)=300.,450.
 mixture_cp_models(1)='polynomial'
 mixture_cp_coefficients(:,1)=3800.,1.,0.,0.
 mixture_cp_temperature_ranges(:,1)=300.,450.
/
```

Coefficients specify `a + b*T + c*T**2 + d*T**3` in SI units with temperature
in kelvin. Defaults are `constant`, using the existing `liquid_density` and
`liquid_cp`. Optional `mixture_vapor_pressure_models` selects `clausius` (default)
or `antoine`; the latter requires `mixture_antoine_coefficients(:,component)`
as `[A,B,C,D]` for `p_sat = D*10**(A-B/(T+C))` and
`mixture_antoine_temperature_ranges(:,component)`. `mixture_boiling_models`
defaults to `reject`; `surface-cap` requires Clausius pressure and a supplied
`mixture_critical_temperatures(component)` above its reference temperature.
This cap is not a physical mixture-boiling rate model.

The shared liquid-property API enforces model/range restrictions. Entries beyond
`component_count` must retain defaults. Omitting the group retains the old
restart context; when provided, the entire group binds checkpoint identity,
including explicitly supplied defaults. Retain the group on restart.
Example fits are software fixtures, not calibrated material data. See
[evidence](../../docs/validation/spray-amr-properties-20260918.md).

### Per-port size distributions

The optional `&amr_spray_sizes` group selects the existing deterministic diameter
sampler independently for each port. For example:

```fortran
&amr_spray_sizes
 port_sizes(1)%distribution='normal',port_sizes(1)%std=2.e-6
 port_sizes(1)%minimum=2.e-6,port_sizes(1)%maximum=3.e-5
 port_sizes(2)%distribution='weibull',port_sizes(2)%shape=2.5
 port_sizes(2)%weighting='number'
 port_sizes(2)%minimum=2.e-6,port_sizes(2)%maximum=4.e-5
/
```

Distributions are `weibull`, `normal`, `lognormal` and `chi-squared`; weighting
is `mass` (default) or `number`. `port_diameters` retains the existing diameter
parameter meaning for the chosen sampler. `std` is the diameter standard
deviation for normal/lognormal fits, in metres; other distributions require
zero. Only Weibull uses a nonzero `shape`. Defaults are Weibull with zero shape,
zero standard deviation and bounds `[1e-12,1] m`, preserving monodisperse
injection. The shared injector validates controls before accepting each event.
Sampling retains deterministic global parcel IDs and event mass accounting.

Omitting the group retains the old restart context. A supplied group binds all
size settings to restart identity; retain it unchanged on continuation. These
example distributions are fixtures, not a calibrated injector. See
[evidence](../../docs/validation/spray-amr-sizes-20260918.md).

### Per-port injection rate histories

`&amr_spray_rates` optionally supplies piecewise-linear mass rates per port:

```fortran
&amr_spray_rates
 profile_points(1)=3
 profile_times(:,1)=0.,1.e-8,2.e-8
 profile_rates(:,1)=.001,.002,.006
 event_end_times=1.e-8,1.5e-8,2.e-8
/
```

Times are seconds, rates kg/s. Each active profile needs 2–64 strictly increasing
knots; supported time/rate bounds are `[0,1e6] s` and `[0,1] kg/s`. For each
event using that port, set its explicit `event_mass` to zero in `&amr_spray`
and provide an end time here. The existing exact piecewise-linear integral
over `[injection_times(event),event_end_times(event)]` becomes its mass.
The interval must lie inside the profile. A positive integral requires a positive
parcel count; a zero integral requires zero parcels. Ports without profiles
retain explicit event masses and require end-time sentinel `-1`.

The entire interval mass is injected at the event's start, using the existing
port speed and parcel count. This is event discretization of a rate history;
divide intervals as needed for temporal resolution. Automatic interval
subdivision is not provided. Optional rate-derived velocity is described below.
Unused profile entries must remain zero and unused end times `-1`. A supplied
group binds the rate table and intervals to restart identity; omitting it
preserves the previous context. See
[evidence](../../docs/validation/spray-amr-rates-20260918.md).

### Rate-derived injection velocity

For ports with rate profiles, an additional group can replace prescribed speed:

```fortran
&amr_spray_velocity
 rate_velocity(1)=.true.
 discharge_coefficients(1)=.8
/
```

The event speed is `mass / (duration * rho_liquid * pi * port_radius**2 * Cd)`.
Density comes from that port's component fractions and injection temperature,
including configured density fits. Set a positive `port_radii` entry in
`&amr_spray`; allowed bounds are `[1e-12,1] m`. The discharge coefficient must
be in `[1e-6,1]` and the resulting speed no greater than `1e7 m/s`. Zero-mass
events use zero speed. The computed scalar speed feeds the existing cone/swirl
geometry. It is an interval-average rate conversion, not a pressure-drop model.

Ports default to prescribed speed and coefficient 1. Non-unit coefficients on
disabled ports and controls on unused ports are rejected. A supplied group
binds restart identity; retain it unchanged on continuation. See
[evidence](../../docs/validation/spray-amr-discharge-20260918.md).

### Rate table files

Each port can load its rate table from a shared text file:

```fortran
&amr_spray_rate_files
 profile_files(1)='injector-rate.dat'
/
```

The file contains exactly two whitespace-separated numeric columns: time in
seconds and mass rate in kg/s. Blank lines and comments starting with `#` or
`!` are allowed. For example:

```text
# time [s] rate [kg/s]
0       0.001
1.e-8   0.002
2.e-8   0.006
```

The same 2–64 knot limits, monotone times and rate bounds apply. Records must
be shorter than 4096 characters. Paths resolve from the process working
directory and must be readable consistently on every rank. Leave that port's
inline `profile_points`, `profile_times` and `profile_rates` at zero; supplying
both a file and inline data is rejected. `event_end_times` still belongs in
`&amr_spray_rates`. The loaded numeric table binds restart identity, not the file
path or comments: moving an unchanged curve does not change the context, while
changing its values rejects restart. Files are read at startup only. See
[evidence](../../docs/validation/spray-rate-files-20260918.md).

### Computational parcel splitting

`options%maximum_multiplicity` in `&amr_spray` now activates the existing
conservative computational splitting after a coupled interval. Zero disables
it. As in the shared spray options, a positive value requires an active
breakup model (`etab` or `khrt`) with its physical inputs. `maximum_parcels`
bounds each owner's resulting parcel count.

Splitting divides represented droplet counts while preserving individual drop
mass, composition and the aggregate inventories. New parcels receive global
IDs and retain parent IDs; the existing checkpoint carries these fields and the
global ID watermark. Capacity failure rejects the entire coupled interval.
Computational child IDs now follow global parent-ID order and sibling order,
independent of owner rank. Future child IDs can differ from older builds that
reserved them in rank order; saved IDs remain unchanged on load.
This does not itself represent additional physical breakup. The coupled API's
optional `split_created` reports the global number of computational children,
excluding KH/splash children created earlier in the interval. See
[focused evidence](../../docs/validation/spray-amr-split-20260918.md).
The application now also has a focused
[2-to-4 rank splitting restart check](../../docs/validation/spray-amr-split-application-20260918.md).

## PPM gas reconstruction

Fixed and selected MPI AMR spray executables accept
`reconstruction='characteristic_ppm'` in the gas namelist. The gas stage uses
four-layer CTU reconstruction with common or recursive level timesteps; the
particle exchange remains in the existing split. Restart context binds the
reconstruction selector. A focused fixed-mechanism periodic 2-to-4-rank
continuation with injection and topology change passed; see the
[PPM application evidence](../../docs/validation/ppm-amr-spray-20260925.md).
Selected runtime and physical spray validation were not part of that check.
