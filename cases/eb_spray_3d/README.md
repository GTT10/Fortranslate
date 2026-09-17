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

`spray.nml` has one `eb_spray` group:

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

Checkpoint/restart settings are explicitly rejected. Independent nozzle ports, physical
outer-boundary particle events, wall-film formation and distributed EB/AMR
particles remain open. Trial trajectories must stay in the root domain. The
particle source/drift split is first order; this synthetic implementation case
is not calibrated spray validation or a convergence study.
