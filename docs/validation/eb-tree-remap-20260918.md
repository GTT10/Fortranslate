# Conservative multilevel EB tree remap — 2026-09-18

Implementation: `acbca12911bd981066d82f07eafc7d9a652448c4`, based on
`2705f28`. API: `remap_reactive_eb_tree_3d`.

## Contract

The source and destination share root cells and refinement ratios on common
levels. The supplied full geometry describes the deeper tree's finest lattice.
The routine coarsens this geometry, privately averages old descendants, and
checks that the reconstructed old-node metrics match the source. It initializes
the new topology and fills cells from the deepest old donor at the same or
coarser level, using integer global coordinates. Same-level overlap copies
state and temperature; new refinement uses ancestor PCM. Old averages supply
coarsening. Final bottom-up averaging precedes publication of the complete
new topology, geometry, state and temperature.

The old tree remains unchanged on topology, geometry, donor or EOS failure.
Moving geometry and changes to common-level refinement ratios are rejected.
Coarsening necessarily loses old subcell detail; subsequent PCM refinement
preserves inventory rather than reconstructing that discarded variation.
The implementation currently builds private trees and full geometry levels;
large-tree memory use and performance are unqualified. Scheduled/adaptive
application integration remains separate work.

## Completed focused verification

GNU Fortran 13.3.0, `/usr/bin/gfortran`, existing Debug/Ninja
`build/selected-tree-mpi-debug`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_eb_tree_remap_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_eb_tree_remap_3d$' --output-on-failure --output-log build/eb-tree-remap-checks/ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, zero failures/skips, 0.58 s (total 0.61 s). Six fixtures
combine x/y/z cuts at 0.31 with four/five-node three-level trees and ratio-2
refinement. Node-dependent states and spatial variation within a grandchild
make overlap preservation nontrivial. Each fixture moves a grandchild,
removes the other grandchild if present, coarsens to root only, then refines
back to the original topology.

Same-level overlapping states and temperatures match exactly. All 18
composite inventory comparisons pass the 3e-14 scaled bound, with maximum
`abs(after-before)/max(1,abs(before))` **7.1245733823814598e-15**. Negative deep
source density and a changed cut position reject with unchanged node count,
parent identities, states and temperatures. No ratio-change, moving-geometry,
application schedule, MPI or convergence qualification is claimed.

Provenance, project contract and whitespace checks pass. No Python source
changed. No failed build/test attempt occurred. Logs:
`build/eb-tree-remap-checks/{build,ctest}.log`; fixture residuals:
`build/selected-tree-mpi-debug/Testing/Temporary/LastTest.log`.
