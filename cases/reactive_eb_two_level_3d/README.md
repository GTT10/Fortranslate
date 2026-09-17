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
Set `automatic_timestep=.true.` to further limit the common interval using
both grids' cut-cell hydro and diffusion estimates. `hydro_cfl` defaults to
0.5 (range `(0,1]`), and `transport_cfl` to 0.35 (range `(0,0.5]`).
Coarse parents are averaged privately before estimation; signal speeds and
diffusivity bounds include donors from both grids. Automatic selection is off
by default, retaining fixed-cap behavior. The selected interval bounds each
transport half-step conservatively. Set `subcycling=.true.` (default false)
to use `refinement_ratio` fine hydro steps and its square in each transport
half. Automatic selection scales fine-grid limits by those step counts while
retaining the coarse-grid limits. Reaction halves operate on composite active
cells around the complete transport/hydro sequence.
[Coupled subcycling evidence](../../docs/validation/eb-subcycle-full-20260917.md). These pre-step estimates do not control
chemistry error or guarantee coupled stability.
[Selection evidence](../../docs/validation/eb-two-level-auto-dt-20260917.md).
Rejected updates halve the trial interval
without changing either grid or accumulated budgets. `maximum_step_retries`
(default 12, range 0–60) limits halvings per accepted step; zero disables retry.
The next step starts again from the user maximum, clipped by remaining time.
Exhaustion stops before CSV output. Rejection-based retry does not establish
CFL stability, temporal error control or physical accuracy.
[Retry evidence](../../docs/validation/eb-two-level-retry-20260917.md).
This application uses plain Euler
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
chemical integrator, `time_step`, automatic-selection mode and CFL controls
and `subcycling` must match the saved run. Common stepping uses V1/V2
checkpoints for fixed/automatic selection; subcycling uses V3/V4 respectively.
These share a record layout with mode-specific magic.

A single file contains both grids and cumulative budgets. The reader validates
the complete file before restoring either grid. A failed read leaves existing
CSV outputs untouched. Checkpoint replacement is staged for validation but is
not atomic against storage failure. Use distinct paths for inputs, outputs,
restart source and checkpoint destination.
[Restart evidence](../../docs/validation/eb-two-level-restart-20260917.md).

The current executable is fixed-mechanism and serial. Selected-mechanism CLI
support, dynamic/multiple patches and MPI are not yet connected.
[Application evidence](../../docs/validation/eb-two-level-application-20260917.md).
