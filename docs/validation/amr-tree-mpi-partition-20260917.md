# MPI cell-step-weighted patch ownership — 2026-09-17

## Source and implementation

Base `7d32500f70416ae87e8c2c12b71d5e2f83a97d21`, local `main`.
Adds deterministic largest-work-first whole-patch assignment with checked int64
weights and transactional owner-vector publication. MPI fixed/selected applications
accept optional `amr_tree_distribution partitioning='weighted'`, defaulting to
`round_robin`. Initialization, remap destinations and restored geometry use the
same policy. Checkpoint identity remains independent of owners and partition policy.
Numerical kernels, mechanisms, tolerances and golden hashes are unchanged.

SHA-256:

- `src/parallel/mpi_amr_patch_tree_distribution_3d_mod.F90`:
  `43db6347b15e9be270cd1358c562f7a5cc543aa76d2fbfb0c036e6e1db2c43ac`
- `src/driver/mpi_amr_patch_tree_application_3d_mod.F90`:
  `51519803cb6357d9f88c23554e132985a13643744e636e51572ff99b997d490f`
- `tests/unit/test_mpi_amr_patch_tree_distribution_3d.F90`:
  `fae2bacae686226133850e75168938ac5841557f79575606f423299909f7642e`
- `tools/check_mpi_amr_tree_adaptation_3d.py`:
  `61d5fb0dfbf4f64518aa41fe29a33041868ce71b63b4d9a577f87744c844c592`

## Commands

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, existing
`build/selected-tree-mpi-debug`; compiler/backend settings unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_mpi_amr_patch_tree_distribution_3d pelef_mpi_amr_tree_reactive_3d test_selected_elementary_mpi_amr_tree_3d -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(unit_mpi_amr_patch_tree_distribution_3d_[23]|regression_mpi_amr_tree_weighted_adaptation_3d)$' --output-on-failure
python3 -m py_compile tools/check_mpi_amr_tree_adaptation_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

## Completed results

Fixed/selected application and distribution test builds passed. CTest **3/3
passed**, 43.18 s total: distribution 2/3 ranks 0.44/0.46 s, weighted adaptation
42.26 s. No failures/skips. Provenance, Python syntax, project contract and
whitespace checks passed. Existing compiler warnings remain.

## Coverage

Distribution unit tests retain prior scatter/gather checks and add deterministic
cell-weight, recursive hydro/transport-weight, tie-breaking, empty-rank and
invalid-rank/overflow rollback assertions. The three-node fixture has cell counts
64,64,27 and ratios 2,3. Recursive transport weights are 64,256,972: weighted
two-rank assignment is [1,1,0], with estimated loads [972,320], versus round-robin
[1036,256]. These are work estimates, not measured speedups. Four-rank recursive
hydro assignment is [2,1,0], leaving one rank empty.

The weighted application test reuses the recursive PLM automatic-adaptation case:
fixed/selected two-rank runs save with weighted owners; three-rank continuation
explicitly switches to round-robin. Serial continuous/restart outputs and checkpoint
bytes must match. The test includes refinement to three levels and ten nodes,
three coarse steps, fixed/selected final byte equality, unchanged output on
rejected mode/threshold/cadence/counter/node-budget inputs, and serial interchange.

Logs: `build/amr-tree-mpi-partition-checks/{build,ctest}.log`.
Artifacts: `build/selected-tree-mpi-debug/tests/mpi_amr_tree_weighted_adaptation_3d/`.

## Limits

Weights estimate cells times cumulative step counts, with transport dominating
when enabled. They do not measure chemistry stiffness, communication, memory or
spatial locality. Whole patches cannot be split. The selection algorithm is
quadratic in patch count plus nodes times ranks; no production scaling is claimed.
Independent runtime load-triggered migration and distributed output remain open.
The application test is hydro-only; recursive transport weighting is checked at
partition API level. No full matrix, installed audit, external solver parity or
physical validation was run.
