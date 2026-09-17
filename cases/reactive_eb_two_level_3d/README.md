# Experimental two-level reacting EB case

```sh
cmake --build build/selected-tree-mpi-debug --target pelef_reactive_eb_two_level_3d
build/selected-tree-mpi-debug/pelef_reactive_eb_two_level_3d cases/reactive_eb_two_level_3d/reference.nml
```

The `eb_two_level` namelist describes a full-domain **finest geometry lattice**
(`nx,ny,nz`). Coarsening by `refinement_ratio` constructs the coarse domain;
`fine_lower` and `fine_upper` are inclusive coarse-cell indices of the one fine
patch. Geometry metrics are extracted from the same finest construction.
The reference uses an 8-cubed lattice, a 4-cubed coarse domain and a 4-cubed
fine patch, with fixed elementary H2/O2 chemistry, molecular transport and hydro.

`geometry_type` accepts `plane`, `extruded` or `level_set`, with the same geometry
file formats as the existing [extruded](../reactive_eb_extruded_3d/README.md) and
[nodal](../reactive_eb_level_set_3d/README.md) cases. Coarsening rejects shapes
that the single-fluid-state coarse geometry cannot represent. `domain_lower`
and `domain_upper` default to the unit cube. Initial densities, pressure and
composition use the existing elementary density-sheet defaults.

`time_step` is a user-supplied fixed maximum interval, clipped at `final_time`.
Choose it to satisfy the explicit small-cell hydro and diffusion limits; there
is no automatic CFL selector or retry loop. This application uses plain Euler
components, with reaction/transport half-steps around hydro. Physical boundaries
are outflow for hydro and zero flux for molecular transport.

`output_file` contains composite active cells only: unrefined coarse cells at
level zero and fine cells at level one. Coordinates are fluid centroids.
`<output_file>.budget.csv` stores initial/final conserved inventories and
accumulated boundary, wall and chemical contributions. The application checks
their balance before writing outputs, using a scale with an absolute floor of one.

## Checkpoint continuation

Set `checkpoint_file` and a positive `checkpoint_interval_steps` to save after
accepted steps. `stop_after_checkpoint=.true.` stops after the first scheduled
save and writes the current CSV outputs. Restart by setting `restart_file`;
`final_time` and `maximum_steps` may increase. Geometry, patch topology, physics,
chemical integrator and `time_step` must match the saved run.

A single file contains both grids and cumulative budgets. The reader validates
the complete file before restoring either grid. A failed read leaves existing
CSV outputs untouched. Checkpoint replacement is staged for validation but is
not atomic against storage failure. Use distinct paths for inputs, outputs,
restart source and checkpoint destination.
[Restart evidence](../../docs/validation/eb-two-level-restart-20260917.md).

The current executable is fixed-mechanism and serial. Selected-mechanism CLI
support, dynamic/multiple patches and MPI are not yet connected.
[Application evidence](../../docs/validation/eb-two-level-application-20260917.md).
