# Experimental serial regular 3D patch-tree application

Build and run from the repository root:

```sh
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_amr_tree_reactive_3d
mkdir -p build/tree-example
cd build/tree-example
../spray-mpi-migration-gnu-debug/pelef_amr_tree_reactive_3d \
  ../../cases/amr_tree_reactive_3d/hotspot.nml \
  ../../cases/amr_tree_reactive_3d/branched.dat
```

The executable accepts `<reactive_3d.nml> <topology.dat> [schedule.dat]`. Flow input uses the
existing `reactive_3d` configuration, including elementary/full-H2/O2 models,
PCM/characteristic-PLM face reconstruction, physical boundaries, chemistry and
molecular transport. Existing validation requires at least four root cells in
each direction. The application is also registered for installation.

## Topology file

The first record is the number of child patches (zero is a root-only run).
Each following record contains:

```text
parent ilower iupper jlower jupper klower kupper ratio
```

Root node is 1. Record `i` defines node `i+1`; its parent must already exist.
Bounds are inclusive parent-cell indices. Siblings cannot overlap. Refinement
ratios must agree at each depth. Counts/depths are dynamically allocated, with
checked coordinate/extent limits. `branched.dat` defines two adjacent root
children and one grandchild: four nodes over three levels.

Each node samples the analytic initial condition at its own global cell centers.
Covered parents are then restricted and temperatures recovered. Time integration
uses synchronous R/2-T/2-H-T/2-R/2, with disabled operators skipped. The timestep
includes all-node hydro/transport bounds, touched isothermal walls and inflow
states at the finest spacing. Rejected steps halve the timestep up to 20
attempts; only successful steps advance time and boundary budgets.

## Scheduled regridding

An optional third argument gives a schedule. The first record is its event count;
each following record is `time 'topology-path'`. Times must increase strictly,
be positive, and not exceed `final_time`. Quote paths, including absolute paths;
relative paths resolve beside the schedule file. All topology inputs are read
before evolution. `regrid.dat` moves/refines the patches with ratio 3, removes
all children, then restores `branched.dat` at final time:

```text
3
2.5e-10 'moved_ratio3.dat'
5.0e-10 'root.dat'
1.0e-9 'branched.dat'
```

Time steps stop exactly at each event. Conserved overlap projection and EOS
recovery occur before the new topology is published. An absolute-cell-inventory
scaled conservation check must pass at `5e-12`; boundary budgets and original
initial inventories are carried through unchanged. Final-time events run before
output. Reservoir bounds use the finest spacing across the initial and scheduled
trees; per-step hydro/diffusion/wall bounds use the active tree. Logs report event
times, node/level counts, and maximum remap error. This is a prescribed schedule,
not automatic solution-based tagging.

## Output

`output_file` is a single CSV of **exposed cells only**, so covered parent cells
are not counted twice. `node,level,dx,dy,dz` precede the standard reactive-3D
columns. Cell volume is `dx*dy*dz`; field integrals must use that weight.
The companion `<output_file>.budget.csv` records `component,initial,final,
boundary_change`. Components are density, x/y/z momentum, total energy, then
species densities in mechanism order. Chemistry changes individual species;
the corresponding source contribution is not a boundary flux.

`hotspot.nml` exercises elementary chemistry, transport and moving no-slip walls.
The small case produces 127 exposed cells and is an implementation regression,
not physical validation. Output write failure may leave a partial CSV; the
files are not checkpoints.

## Remaining scope

This frontend supports static or scheduled topologies with serial ownership.
Automatic tagging, checkpoint/restart, MPI execution and a selected-mechanism frontend
remain to be connected. Coarse/fine donors are PCM and time stepping is common
across levels; this is not subcycled or globally second-order AMR qualification.
