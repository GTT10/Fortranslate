# Synchronized two-level EB hydro — 2026-09-17

Implementation `8a3b409a04ffe698c27262ab6f23240d143e1dcb`, based on
`911c7c0` on local `main`.

## Implemented operation

`advance_reactive_eb_two_level_euler_3d` advances one static fine patch and its
coarse domain through a common Euler interval with outflow physical boundaries.
The patch may be internal or touch the physical domain. It verifies the fine
geometry restricted to coarse resolution against the corresponding coarse box,
including classification, volume, faces, centroids, wall area and normals.

The operation privately averages initial fine states into refined coarse cells,
constructs PCM coarse-donor face halos, advances both grids, restricts the accepted
fine flux, replaces coarse patch-face fluxes, refluxes adjacent unrefined coarse
cells and averages the accepted fine state down. Existing EOS validation is used
at each state-changing stage. Coarse/fine conserved states, temperatures and both
returned inventories are committed together only after all stages succeed.

Returned `boundary_change` is the per-step composite outer-face flux integral;
faces covered by the patch use fine-restricted flux. `wall_change` is the per-step
EB wall source integrated over unrefined coarse cells and all fine cells. Their
sum accounts for the composite state change. They overwrite the supplied output
vectors on success and remain unchanged on failure.

This is a serial library integrator with one static patch and hydro only.
It has no chemistry/transport splitting, subcycling, time interpolation, dynamic
regridding, multiple patches, hierarchy checkpoint/application or MPI connection.
It uses plain Euler; callers must supply a suitable small-cell-limited timestep.
Nonphysical updates are rejected, not stabilized or clipped. Existing single-level
application behavior is unchanged. This increment does not establish full PeleC
parity or completion of EB-AMR support.

## Independently completed evidence

GNU Fortran 13.3.0 Debug, Open MPI 4.1.6, Ninja; existing matching build
`build/selected-tree-mpi-debug`, MPI enabled, spray disabled. Core and selected
EB runtime compiled. Execution uses fixed elementary thermodynamics; no MPI or
application matrix was run. No Python source changed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 2 --target test_eb_two_level_3d pelef_selected_reactive_eb_3d_runtime
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 2 --target test_eb_two_level_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_eb_two_level_3d$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **1/1 PASS, 0 failures, 0 skips**, test 0.30 s, total 0.33 s.
A 4-cubed coarse domain with a cut plane uses ratios two and three, an internal
patch and a patch touching physical boundaries. Every configuration completes
three nontrivial Rusanov steps at 1e-10 s. Each conserved component, including all
species, satisfies the composite boundary-plus-wall inventory balance. Refined
parent cells match child inventories cell by cell, temperatures stay positive,
and covered coarse storage is unchanged. Excessive timestep, inconsistent fine
geometry and invalid ratio reject with all four state arrays and both inventories
unchanged. This is an actual two-level update, not a synthetic reflux-only test.

Provenance and source contract passed. Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Logs: `build/eb-two-level-checks/{build,final-build,ctest}.log`.
An earlier pass preceded the explicit parent/child inventory assertion; only the
final run is counted. Convergence, curved-geometry hierarchy runs, external solver
parity and physical validation are not established by this focused check.
