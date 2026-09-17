# Temperature-adaptive multilevel EB application — 2026-09-18

Implementation: `162377145ab0869f47ed98889b9ff1e8d73c29db`, based on
`c172fff`. [Case guide](../../cases/reactive_eb_tree_3d/README.md).

The application now evaluates exposed-fluid temperature tagging at prescribed
evaluation times and conservatively remaps before saving checkpoints. Inputs
specify requested depth (1:8), per-level ratios/thresholds, tile shape, buffer
and node capacity (1:33, including root). Prescribed event topology is rejected
in adaptive mode. Geometry is prebuilt on the deepest initial/requested lattice.

V4 checkpoints bind tagging policy, evaluation times and initial topology,
store completed generated layouts, and require absent future layouts. Readers
validate topology reconstruction, capacity, depth/ratios, parent ordering,
tile alignment/extents/order, embedded physics/geometry/state and stage/time.
Past temperature decisions are not recomputed. Initial topology may differ
from the generated tiling. This does not add selected frontend, subcycling,
MPI, large-tree scaling or physical qualification.

## Completed focused verification

GNU Fortran 13.3.0, `/usr/bin/gfortran`, existing Debug/Ninja
`build/selected-tree-mpi-debug`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/check_eb_tree_application_3d.py
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target pelef_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_eb_tree_(scheduled|adaptive)_3d$' --output-on-failure --output-log build/eb-tree-adaptive-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_eb_tree_adaptive_3d$' --output-on-failure --output-log build/eb-tree-adaptive-checks/ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

Two distinct regressions passed: scheduled **1/1, 6.41 s** in the initial
combined invocation; adaptive final **1/1, 10.81 s** (total 10.84 s) after a
test expectation correction. No skips. No fresh combined all-pass invocation
or broader matrix is claimed.

The initial five-node three-level state ranges from approximately 365 to
1095 K. Thresholds 1000/1e9 select four first-level children at both events,
leaving cold root cells exposed: **272 exposed cells** at levels 0/1. The
maximum normalized inventory residual is **8.32667281161794579e-16**;
independent CSV volume sums agree with final budgets. Restart before event
one and after either event yields byte-identical final state and budget CSVs.
Thresholds 1e9/1e9 coarsen to root only. Capacity one, altered tag buffer,
changed evaluation time, inconsistent stage and malformed records reject.
Existing protected-file checks also pass.

The first adaptive assertion incorrectly expected only level-1 output, assuming
all fluid exceeded 1000 K. Live CSV temperatures showed 144/272 exposed cells
above the threshold before adaptation; the expected level set was corrected
to {0,1}. Solver behavior and numerical tolerances were unchanged. That failure
and the scheduled pass are retained in `initial-ctest.log`.

Provenance, Python syntax, project contract and whitespace checks pass. Logs:
`build/eb-tree-adaptive-checks/{build,initial-ctest,ctest}.log`. Artifacts:
`build/selected-tree-mpi-debug/tests/eb_tree_adaptive_3d/`. The application
fixture exercises refinement removal/first-level selection, not creation of
new deepest-level patches; the preceding tagging unit fixture covers fine-only
peaks selecting two refinement levels.
