# Periodic MPI mixed spray

Build with `PELEF_ENABLE_MPI=ON` and `PELEF_ENABLE_SPRAY=ON`, then run from a
fresh output directory (input paths may be absolute):

```sh
mpiexec -n 4 /path/to/build/pelef_mpi_spray_3d /path/to/cases/mpi_spray_3d/periodic_mixture.nml
```

The fixed H2/O2 and selected-mechanism applications advance periodic x-slab gas and two-component
particles with P/2-C/2-T/2-H-T/2-C/2-P/2 splitting. The example uses synthetic
water/peroxide liquid properties; it is an implementation exercise, not a qualified
physical spray model. Global cells must provide at least two x cells per rank and
three y/z cells. Rank zero reads and broadcasts the `mpi_spray` namelist.

The default initial gas is 800 K, density 1 kg/m3 and mass fractions H2=0.03,
O2=0.23, N2=0.74. `gas_fractions` follows the loaded mechanism order. The fixed order is
H2,H,O,O2,OH,H2O,HO2,H2O2,AR,N2. Default fractions are assigned by the names H2,
O2 and N2, so their positions may change in a selected mechanism. Liquid species, densities, heat capacities,
boiling temperatures (at 101325 Pa), latent heats and injection mass fractions
are configurable pairs. `number=0,total_mass=0` disables liquid injection.
`injection_events` (default 1), `injection_start` (default 0) and
`injection_interval` (default 1e-8 s) define repeated cone events. `number` and
`total_mass` apply to each event. Setting `injection_events=0` disables all events.
Steps end exactly at pending event times; events at the final time are processed
before checkpoint/output. Each event preserves prior particles and reserves new
IDs, then migrates the new population to owners.

Time controls are `final_time`, `dt_max`, `cfl`, `max_steps`, `rtol`, `atol` and
`chemistry` (`implicit`, `explicit`, or `cvode` in a SUNDIALS build). The code
reduces hydro/transport time-step limits across ranks and retries rejected complete
updates with halved steps. `transport_flags(4)` selects viscosity, thermal
conduction, species diffusion and barodiffusion. Solver/limiter defaults are HLLC/MC.

Outputs use `output_prefix`:

- `.history.csv`: step, time, global parcel count, particle substeps, total
  gas-plus-liquid mass/momentum/energy, relative mass/energy drift after subtracting
  accumulated injection, event count, current ID watermark and injected budgets. Comment
  headers retain the mechanism SHA and initial ID high-water mark.
- `.rankNNNNNN.gas.csv`: local cells with global integer indices, temperature and
  conserved gas fields.
- `.rankNNNNNN.particles.csv`: IDs, position, velocity, per-drop mass/temperature,
  multiplicity and per-drop component masses in configured liquid order.

Existing files are not overwritten. CSV output is not a restart checkpoint.
Distributed KH child release, wall films and SGS remain to be integrated in this executable.

## Save and resume

Set `checkpoint_path='step.chk'` to save the final accepted state. Set
`restart_path='step.chk'` in a subsequent run, choose a fresh `output_prefix`, and
increase `final_time` to continue. Already completed injection events are skipped on restart. Gas,
particles, component masses, budgets, time, step/event counts and the ID high-water mark
are restored. The MPI rank count may change if the new partition still has at
least two x cells per rank. The implementation gathers checkpoint data on root;
it is intended for bounded cases, not production parallel I/O.

The checkpoint binds the mechanism, mixture evolution, geometry, liquid properties,
injection definition and numerical controls. These must match. Rank count,
checkpoint/output paths, `final_time` and `max_steps` may change. Invalid or
incompatible checkpoints are rejected before output history is created. Existing
checkpoint destinations are refused. Focused same-rank continuation matches final
CSV values exactly; changed-rank continuation is compared with floating-point
tolerance, not claimed bitwise identical.


The scheduled-injection context is `mpi-spray-application-v2`: reference budgets
contain the initial gas and accumulated injection is stored separately. Older v1
application checkpoints are rejected because they used a different reference
budget convention. The underlying checkpoint file schema is unchanged.

## Selected mechanism

When `PELEF_MECHANISM_BUNDLE` is set alongside MPI and spray, CMake also builds
`pelef_mpi_spray_3d_selected`. It uses the same application and checkpoint code
with mechanism-sized gas arrays and the selected bundle SHA in checkpoint context.
The existing importer/provenance and selected thermodynamic/transport validation
remain in effect. Liquid vapor names must exist in the selected mechanism.

`periodic_methanol_water.nml` exercises the 21-species FFCM1_Red bundle prepared
by the [spray mechanism recipe](../../docs/spray.md). Its water/methanol liquid
properties are synthetic defaults. Example:

```sh
mpiexec -n 4 /path/to/build/pelef_mpi_spray_3d_selected /path/to/cases/mpi_spray_3d/periodic_methanol_water.nml
```

The old selected MPI 1D executable requires its own 2/10-species initialization
profiles and is omitted for other species counts. This does not disable the
selected MPI spray or AMR targets. Qualification here covers 21 species, not every
supported mechanism size or physical liquid-fuel ignition.

## Time-varying mass flow

Set `injection_profile_points` (2--64), `injection_profile_times` (s) and
`injection_profile_rates` (kg/s, 0--1). Times must increase and cover the complete
window from `injection_start` through
`injection_start + injection_events * injection_interval`; unused array entries
must stay zero. Each event uses the exact piecewise-linear integral over its
interval as the injected mass, replacing `total_mass`. The mass for that interval
is introduced at its event's start, consistent with the discrete parcel schedule.
The final event interval is included even if the simulation stops earlier.

A zero integrated dose consumes the event but creates no particle or ID. Positive
doses use `number` particles. Profiles and their integration version are bound to
checkpoint context, so a changed profile cannot silently continue a saved run.
With `injection_profile_points=0` (default), fixed mass per event is unchanged.


## Independent injection ports

Set `port_count=1..16` to replace the scalar schedule with independent ports.
For each active port, set `port_events`, `port_start` (s), `port_interval` (s),
`port_number` and `port_mass` (kg/event). Defaults are zero events, start zero,
interval 1e-8 s, four parcels and 1e-10 kg. The combined event count is at most
one million. Events execute in increasing time; ties execute in port-index order.
A checkpoint taken at an event time includes every event due at that time.

Port properties are `port_center(:,p)`, `port_axis(:,p)`, `port_speed(p)`,
`port_half_angle(p)`, `port_diameter(p)`, `port_temperature(p)` and
`port_fractions(:,p)`. Their factory defaults match the scalar defaults; setting
scalar inputs does not override port defaults. The two liquid component property
tables remain shared. For example:

```fortran
 port_count=2
 port_events=2,2
 port_start=0,7.5e-9
 port_interval=7.5e-9,7.5e-9
 port_number=2,3
 port_mass=1e-10,2e-10
 port_center(:,1)=.002,.002,.002
 port_center(:,2)=.006,.002,.002
 port_axis(:,2)=-1,0,0
 port_fractions(:,1)=1,0
 port_fractions(:,2)=0,1
```

An enabled shared rate profile applies separately to each port's event interval
and overrides `port_mass`. It must cover every port interval. Overlapping intervals
inject separate doses, so the supplied rate is per port, not a combined flow rate.
Each port can override this shared profile with its own table (below).
Active port definitions are bound to checkpoint context. Rank-count changes retain
completed event order and ID reservation. `port_count=0` preserves scalar behavior
and the existing scalar checkpoint context.


### Port-specific flow histories

`port_profile_points(p)=2..64` enables `port_profile_times(:,p)` (s) and
`port_profile_rates(:,p)` (kg/s, 0--1) for port `p`. Each table must increase in
time and cover that port's complete event intervals. Trailing unused entries
must be zero. Profiles for inactive ports or ports with no events are rejected.

A port-specific table overrides the shared profile for that port. Ports with
`port_profile_points(p)=0` use the shared profile when enabled, otherwise their
fixed `port_mass`. The shared table need only cover the ports that use it. For
example, added to the two-port example above:

```fortran
 injection_profile_points=2
 injection_profile_times=0,1.5e-8
 injection_profile_rates=.01,.01
 port_profile_points=0,3
 port_profile_times(:,2)=7.5e-9,1.5e-8,2.25e-8
 port_profile_rates(:,2)=0,0,.04
```

Port 1 receives two 7.5e-11 kg doses. Port 2 consumes its first event without
particles, then receives 1.5e-10 kg. The checkpoint binds every supplied table.
The checkpoint context capacity is now 131072 characters to retain all 16 tables
at 64 points each; existing shorter contexts and file schemas are unchanged.


## Nozzle geometry and size distributions

The MPI application now passes the existing deterministic injector options through
to scalar or independent-port injection. Namelist derived-type fields select the
same models as the serial injector:

```fortran
 nozzle_radius=4e-4
 nozzle%profile='solid'
 nozzle%swirl=0.1
 size_options%distribution='lognormal'
 size_options%std=4e-6
 size_options%minimum=1e-5
 size_options%maximum=4e-5
 size_options%weighting='mass'
```

With independent ports use `port_nozzle_radius(p)`, `port_nozzle(p)%profile`,
`port_nozzle(p)%swirl`, `port_nozzle(p)%spread` and `port_size(p)%...` instead.
Port defaults are independent of scalar settings.

Nozzle profiles are `legacy` (default), `solid` and `hollow`. The latter two
require a positive outlet radius (m). `swirl` and the hollow-only `spread` are
angles in radians. The outlet lies normal to the port axis at its center.

Size distributions are `weibull` (default), `normal`, `lognormal` and
`chi-squared`, with `mass` or `number` CDF weighting. `minimum`/`maximum` are
truncation diameters (m), default 1e-12/1. `shape` is the Weibull exponent (default
zero gives a single size); `std` is used by normal/lognormal distributions in m.
Defaults preserve the previous monodisperse point-cone behavior. The requested
`diameter` / `port_diameter` provides the size parameter; see the serial
[spray model guide](../../docs/spray.md) for distribution definitions.

Each nondefault injector configuration enters checkpoint identity. Changing it
rejects restart. These are parcel-source models, not resolved nozzle geometry or
experimental injector qualification.


## Representative-drop breakup

Set `breakup_model='etab'` or `'khrt'`, `surface_tension` (N/m) and
`liquid_viscosity` (Pa s) to enable the existing representative breakup models.
The default is `'none'`. KH-RT constants are `khrt_b0=.61`, `khrt_b1=7` and
`khrt_c3=1`. The same settings apply to every port. For example:

```fortran
 breakup_model='khrt'
 surface_tension=.072
 liquid_viscosity=.001
```

Breakup decreases physical droplet mass and increases multiplicity, preserving
parcel component inventories. Mixture density and thermal properties follow the
current component masses and temperature; surface tension and viscosity are the
supplied common constants. KH child shedding is not enabled in this executable.

Particle CSVs include `deformation`, `deformation_rate` and `rt_time`.
Checkpoints retain the appropriate deformation/RT state together with component
masses and IDs. All breakup parameters are bound to restart identity. Changed
parameters are rejected; default-disabled runs retain their earlier checkpoint
context and schema.
