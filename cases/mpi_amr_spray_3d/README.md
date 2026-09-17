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

For timed injection, set `injection_events` (1 through 1024), strictly increasing
nonnegative `injection_times`, and matching `event_parcels` and `event_mass`
arrays. Each event uses the same liquid composition and cone settings. This
explicit schedule replaces the default time-zero injection; when
`injection_events=0` (default), `parcel_count` and `total_mass` retain their
original time-zero meaning. A zero-mass event must also have zero parcels.
For example:

```fortran
 injection_events=2, injection_times=0.,1.e-8
 event_parcels=4,4, event_mass=4.e-11,6.e-11
```

The solver lands exactly on each event. If regridding and injection coincide,
the new grid is established first. Events at final time are applied before
output; later events remain unapplied. New parcels, component masses, IDs and
the cumulative injection ledger publish only after injection and migration
succeed. The summary reports the number of applied events.

An optional single scheduled regrid is specified by `regrid_time` and
`regrid_topology` in `amr_spray`. The latter path is relative to the run directory
unless absolute. The solver lands on that time and atomically remaps gas and
particles. The root physical domain/resolution remains fixed. `root.dat` is a
coarsening target; `tree.dat` refines one corner.

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
boundaries, initial or scheduled single-cone injection and one optional prescribed
regrid. Multiport injection, automatic tagging/rebalancing, checkpoint/restart, wall
films and EB handling are not connected here. Gas-only AMR checkpoint/adaptation
namelists are not controls for this executable. Use the `amr_spray` controls above.
