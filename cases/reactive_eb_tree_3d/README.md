# Serial multilevel EB reacting tree

Build and run the fixed H2/O2 application:

```sh
cmake --build build/debug --target pelef_reactive_eb_tree_3d
./build/debug/pelef_reactive_eb_tree_3d cases/reactive_eb_tree_3d/input.nml
```

The supplied case has a 4-cubed root, two children, and grandchildren with
separate parents (three levels, ratios two). A planar wall cuts z=0.31.
Two 1e-9 steps run reaction, molecular transport and HLLC hydro.

## Input

The `eb_tree` namelist uses `root_shape(3)` and `child_count` (0:32).
Children are listed in parent-before-child order: child entry p is tree node
p+1, root is node 1. `parent(p)` names its parent node. `lower(:,p)` and
`upper(:,p)` are inclusive cell indices in that parent; `ratio(p)` defines
refinement. Siblings cannot overlap and each level uses one refinement ratio.
An empty child list runs a root-only tree.

`domain_lower/upper` set the physical domain. `geometry_type`, `plane_axis`,
`plane_position`, `extrusion_axis`, `section_file` and `level_set_file` feed
the existing geometry builder on the full finest lattice. Coarsening must
represent the geometry at every level; unsupported geometry is rejected.

`time_step` is the maximum common interval; `final_time` clips the last step.
With `automatic_timestep=.true.`, `hydro_cfl` (default 0.5) and
`transport_cfl` (default 0.35) bound the common interval using every patch's
cut-cell volume/open-face metrics and global donor signal/diffusivity maxima.
The estimate uses an averaged snapshot, includes hidden parents which are
advanced before restriction, and remains capped by `time_step`/remaining time.
The full-interval diffusion bound is conservatively applied despite split
half steps. It is not a chemistry accuracy estimate.
Rejected steps halve their interval up to `maximum_step_retries` (default 12).
`maximum_steps` bounds accepted steps; reaching it early is an error.
`chemistry_enabled` and `transport_enabled` are independent. Enabled transport
uses viscosity, conduction, species diffusion and barodiffusion.
`chemistry_integrator` accepts explicit/implicit, with `rtol/atol` for reaction.
`riemann_solver` accepts rusanov/hllc/pelec. Initial conditions use the existing
fixed-mechanism density-sheet initializer and its defaults.

## Scheduled topology changes

`scheduled.nml` demonstrates 5 nodes -> root only -> 4 nodes at times
1.5e-9 and 2.5e-9. `regrid_count` (0:16), increasing positive `regrid_times`,
`regrid_child_count(event)`, `regrid_parent/ratio(:,event)` and
`regrid_lower/upper(:,:,event)` describe complete replacement topologies.
Zero children coarsens to the root. All layouts share root cells and refinement
ratios on common levels. Geometry is built once on the deepest scheduled lattice.

Steps land on event times before conservative remapping and checkpoint saving.
A V3 checkpoint binds the original topology, every scheduled time/topology and
completed stage. Restart reconstructs that stage's geometry before loading
node records, then checks stage/time consistency. Future schedule edits reject
restart. Events beyond the current final time may remain pending.

## Temperature adaptation

`adaptive.nml` evaluates exposed fluid temperatures at the scheduled
`regrid_times`. Set `temperature_adaptive=.true.` and leave prescribed event
boxes/counts at their defaults. `tag_depth` (1:8), `tag_ratios` and
nondecreasing positive `tag_thresholds` define requested refinement levels;
`tag_tile_shape`, `tag_buffer` and `tag_maximum_nodes` (1:33, root included)
control deterministic tiling. The capacity applies before solid-only tiles
are pruned. Exceeding it rejects the adaptation; hot regions are not dropped.

Geometry covers the deepest initial/requested level. Fine hot peaks vote
directly; child-covered parents and solid temperatures do not. Cold states
may coarsen to the root. V4 checkpoints store the complete tagging policy,
evaluation times and all completed generated layouts. Future layouts remain
absent. Restart validates node bounds, parent order, refinement ratios, tile
alignment/extents/order and saved stage/time before continuing. Completed
layouts are restored, without rerunning past temperature decisions.

## Output and restart

`output_file` contains exposed fluid cells only: node, level, local indices,
physical fluid volume, temperature and conserved components. The adjacent
`.budget.csv` gives initial/final composite inventories, boundary, wall,
chemical exchange and residual for each component. Fluid-volume-weighted
cell sums recover the final inventories without counting covered parents.

Set `checkpoint_file`, `checkpoint_interval_steps`, and optionally
`stop_after_checkpoint=.true.`. Continue with `restart_file` and the same
topology, geometry, physics, integrator, time_step, automatic mode and CFL coefficients. Output names and final
limits may change within the existing embedded checkpoint contract. Use a
different checkpoint destination from the restart source.

The application checks the policy envelope, all node records, topology,
geometry, EOS, parent averages, composite budget and EOF before producing
output. Saving serializes to a scratch stream before opening the destination.
Replacement is not filesystem-atomic if an I/O error occurs during copying.

## Scope

This is a serial static, scheduled or temperature-adapted, fixed-mechanism application using common-interval
Euler operators, spatial PCM donors, hydro outflow and closed diffusion
boundaries. It is not a convergence or physical validation case. Selected
mechanisms, recursive subcycling, MPI and
large-tree scaling remain open for this EB tree application. Plain Euler
retains the small-cell interval restriction; retries do not replace an
accuracy/convergence study.
