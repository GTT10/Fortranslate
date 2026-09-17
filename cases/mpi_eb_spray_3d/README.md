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
injection. Checkpoint/restart, regridding, rebalancing and child creation are
explicitly rejected. It inherits root-contained specular EB drift and the
composite stencil approximation. It is a synthetic integration fixture, not
nozzle, diesel ignition, convergence or physical validation.
