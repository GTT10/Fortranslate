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
periodic domain, chemistry/transport, numerical scheme, CFL and final time.
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
recursive gas steps; the driver retains a conservative finest-grid CFL bound.

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
It is not physical spray validation. This initial executable supports periodic
boundaries, initial or scheduled multiport injection and one optional prescribed
regrid or temperature-driven adaptation, and gas-work-based rebalancing. Wall
films and EB handling are not connected here. Gas-only AMR checkpoint/adaptation
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
