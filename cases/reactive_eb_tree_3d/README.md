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
Rejected steps halve their interval up to `maximum_step_retries` (default 12).
`maximum_steps` bounds accepted steps; reaching it early is an error.
`chemistry_enabled` and `transport_enabled` are independent. Enabled transport
uses viscosity, conduction, species diffusion and barodiffusion.
`chemistry_integrator` accepts explicit/implicit, with `rtol/atol` for reaction.
`riemann_solver` accepts rusanov/hllc/pelec. Initial conditions use the existing
fixed-mechanism density-sheet initializer and its defaults.

## Output and restart

`output_file` contains exposed fluid cells only: node, level, local indices,
physical fluid volume, temperature and conserved components. The adjacent
`.budget.csv` gives initial/final composite inventories, boundary, wall,
chemical exchange and residual for each component. Fluid-volume-weighted
cell sums recover the final inventories without counting covered parents.

Set `checkpoint_file`, `checkpoint_interval_steps`, and optionally
`stop_after_checkpoint=.true.`. Continue with `restart_file` and the same
topology, geometry, physics, integrator and time_step. Output names and final
limits may change within the existing embedded checkpoint contract. Use a
different checkpoint destination from the restart source.

The application checks the policy envelope, all node records, topology,
geometry, EOS, parent averages, composite budget and EOF before producing
output. Saving serializes to a scratch stream before opening the destination.
Replacement is not filesystem-atomic if an I/O error occurs during copying.

## Scope

This is a serial static-tree, fixed-mechanism application using common-interval
Euler operators, spatial PCM donors, hydro outflow and closed diffusion
boundaries. It is not a convergence or physical validation case. Selected
mechanisms, automatic intervals, regridding, recursive subcycling, MPI and
large-tree scaling remain open for this EB tree application. Plain Euler
retains the small-cell interval restriction; retries do not replace an
accuracy/convergence study.
