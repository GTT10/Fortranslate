# Recursive multilevel EB hydro — 2026-09-18

Implementation: `bd0e2817bfbd4d1bedc903d5b05511d6cccfc46c`, based on
`ecadac6`. API: `advance_reactive_eb_tree_subcycle_hydro_3d`.

Each level advances its nodes together with Euler hydro, reconciles shared
peer transfers, and supplies provisional old/new states for time-interpolated
ancestor halos. The next level takes refinement-ratio substeps. Accepted
effective fine transfers accumulate over the parent interval; interval-local
coarse/fine reflux and averaging run before the next enclosing substep.
Ancestor clocks retain independent interpolation fractions. Last fine steps
use the remaining interval to avoid duration-sum drift.

Flux synchronization now optionally restricts its correction/propagation to
an inclusive level range. The recursive algorithm synchronizes only a level
and its immediate children, avoiding reapplication of already completed
deeper corrections. Default whole-tree synchronization is unchanged.

The outer transaction preserves the caller's tree, boundary vector and wall
vector on failure. Wall contributions mask child-covered cells on each level
step; final physical exchange comes from the synchronized root flux. This is
a serial hydro library path. Subcycled molecular transport, reaction/transport
composition, automatic coarse-interval selection, application/restart mode,
MPI, convergence and physical qualification remain separate work. Repeated
private tree copies and halo averaging have not been optimized for scale.

## Completed focused verification

GNU Fortran 13.3.0, `/usr/bin/gfortran`, existing Debug/Ninja
`build/selected-tree-mpi-debug`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_eb_tree_subcycle_hydro_3d test_eb_tree_flux_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_eb_tree_(subcycle_hydro|flux)_3d$' --output-on-failure --output-log build/eb-tree-subcycle-hydro-checks/ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, zero failures/skips, total 1.61 s (flux 0.15 s, recursive
hydro 1.42 s). Hydro fixtures combine x/y/z planar cuts at 0.31 with four/five
nodes, two ratio-2 refinement levels, one-sided deep refinement and peers with
different parents. Nonuniform states advance over a 1e-8 root interval with
HLLC. Composite change minus physical boundary and exposed wall transfers has
maximum normalized residual **6.5088007743307606e-15**, below the unchanged
3e-14 test bound. Deep states change. Invalid solver and invalid deep density
preserve all states, temperatures and sentinel exchange vectors. The existing
whole-tree manufactured-flux test also passes after the range extension.

Provenance, project contract and whitespace checks pass. No Python source
changed. No failed build/test attempt occurred. Logs:
`build/eb-tree-subcycle-hydro-checks/{build,ctest}.log`. Per-fixture output:
`build/selected-tree-mpi-debug/Testing/Temporary/LastTest.log`.
