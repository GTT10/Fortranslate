# Experimental two-level reacting EB case

```sh
cmake --build build/selected-tree-mpi-debug --target pelef_reactive_eb_two_level_3d
build/selected-tree-mpi-debug/pelef_reactive_eb_two_level_3d cases/reactive_eb_two_level_3d/reference.nml
```

The `eb_two_level` namelist describes a full-domain **finest geometry lattice**
(`nx,ny,nz`). Coarsening by `refinement_ratio` constructs the coarse domain;
`fine_lower` and `fine_upper` are inclusive coarse-cell indices of the fine
patch when `patch_count=1` (default). Geometry metrics are extracted from the same finest construction.
The reference uses an 8-cubed lattice, a 4-cubed coarse domain and a 4-cubed
fine patch, with fixed elementary H2/O2 chemistry, molecular transport and hydro.

For multiple static, nonoverlapping fine patches, use
[`multipatch.nml`](multipatch.nml). Set `patch_count` (2–16) and provide
`patch_lower(:,p)` / `patch_upper(:,p)` as inclusive coarse-cell boxes.
Namelist values run through the three coordinates of each patch in order.
These boxes replace `fine_lower` / `fine_upper` in multiple-patch mode.
All patches share the refinement ratio. `subcycling=.true.` advances every
fine patch together through `refinement_ratio` hydro substeps and its square
in each transport half. Peer halos use the same fine-time snapshot, shared
fluxes reconcile after every substep, and accumulated fluxes drive the final
coarse reflux. Automatic selection includes every patch's bounds and donors
and scales fine bounds by the corresponding substep count.
[Multiple-patch subcycling evidence](../../docs/validation/eb-multipatch-subcycle-20260917.md).

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
With multiple patches, an extra leading `patch` column identifies coarse cells
as 0 and fine patches by input index (1 through `patch_count`). Cell indices
remain local to each patch. Budgets include the complete composite domain.
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
Multiple-patch checkpoints use separate `PELEF_EB_MULTIPATCH_CHECKPOINT_V1`
(fixed) or `V2` (automatic) headers for common stepping; subcycling uses `V3`
or `V4` respectively. They include all fine records and require the same patch
order, boxes and stepping mode on restart.

A single file contains all grids and cumulative budgets. The reader validates
the complete file before restoring any grid. A failed read leaves existing
CSV outputs untouched. Checkpoint replacement is staged for validation but is
not atomic against storage failure. Use distinct paths for inputs, outputs,
restart source and checkpoint destination.
[Restart evidence](../../docs/validation/eb-two-level-restart-20260917.md).

Both executables are serial. Dynamic patches and MPI are not yet connected.
[Multiple-patch application evidence](../../docs/validation/eb-multipatch-application-20260917.md).
[Application evidence](../../docs/validation/eb-two-level-application-20260917.md).

## Selected mechanism

Configuring `PELEF_MECHANISM_BUNDLE` also builds and installs
`pelef_reactive_eb_two_level_3d_selected`. It shares the application source and
supports the same single/multiple-patch paths. Supply `composition_count`,
`composition_species` and `composition_mole_fractions`; species names resolve
against the pinned bundle and fractions are normalized by the existing resolver.
The default chemistry integrator comes from the bundle; explicit and implicit
integrators can be selected in the namelist. Fixed-mechanism runs reject
composition overrides.

With `mechanisms/h2o2_full.json`, run
[`selected_h2o2.nml`](selected_h2o2.nml). This example uses two fine patches,
automatic interval selection, ten-species chemistry and molecular transport.
Checkpoint patch records bind bundle SHA, resolved composition and integrator,
so selected/fixed or changed-mechanism/composition restarts are rejected.
The initial regular/cut temperatures must lie within the common species
thermodynamic range. The inherited elemental inventory supports only
H2, H, O, O2, OH, H2O, HO2, H2O2, N2 and AR; other species are rejected by this
application, including when reaction is disabled. This is not a general
carbon-fuel EB frontend or a CVODE connection.
[Selected application evidence](../../docs/validation/eb-selected-two-level-20260917.md).
