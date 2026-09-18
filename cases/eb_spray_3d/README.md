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

Independent nozzle ports, physical outer-boundary particle events and distributed
EB impact/film ownership remain open. Trial trajectories must stay in the root domain. The
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
reactive update. KH shed-mass child creation remains unsupported in this serial
application. See [focused evidence](../../docs/validation/spray-eb-split-20260918.md).
