# Selected-level MPI EB SDC loop — 2026-09-27

Implementation `f5ae4979`, based on `4fc25017` on local `main`.
PeleC pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`advance_mpi_eb_tree_sdc_3d` now accepts zero-based level selection, ancestor
state/temperature endpoints, start/end gas fractions and predictor-source
endpoints. Source halos use the midpoint of the supplied start/end fractions.
MPI fractions remain root-first and sized to the maximum tree level. Optional
control groups, level choice and requested outputs are agreed collectively.

Selected-level mode skips initial restriction and final reflux. It advances
only that level, retains other gas/temperature/reaction estimates, and requires
an integrated-flux output for later synchronization. Source/wall/chemical
inventories include exposed cells of the selected level only; physical-boundary
inventory is returned only for root-level or common stepping. Optional old/new
predictor outputs include external, transport and reaction rates, with inactive
levels zero. All outputs publish only after collective success. Common-step
behavior remains available through omitted level selection.

## Completed focused checks

GNU Fortran 13.3.0, `/usr/bin/mpifort` and `/usr/bin/mpiexec`
(OpenMPI/OpenRTE 4.1.6), Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance and the touched test target built successfully.
- CTest `unit_mpi_eb_tree_full_3d_2`: 1/1 passed,
  78.93 s test / 78.98 s total.
- New comparisons use three levels, five patches, depth 2, dt 1e-10 s,
  1000 K initial fluid, two SDC iterations, chemistry, all transport flags,
  external forcing and outflow boundaries. Both round-robin ownership and
  all work on rank 0 (rank 1 idle) are tested.
- Gas end states are 1.01 times their start values; start/end fractions are
  0.25/0.75. Predictor source endpoints also differ. Selected state, temperature,
  reaction estimates, all face-flux arrays and old/new predictor fields match
  serial exactly. Limiter is exact; inventories satisfy the asserted 2e-12
  relative agreement bound.
- Inactive gas/temperature remain unchanged, inactive reaction sentinels remain
  7, predictor rates are zero and inactive flux arrays remain unallocated.
- Invalid artificial viscosity rejects after old transport without publishing
  state, temperature, reaction estimates, inventories, fluxes or predictors.
- Existing common-step, reaction-stage and split-step checks remain in the
  executable and passed. Python syntax, project contract and diff checks passed.

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_full_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_full_3d_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/sdc-mpi-eb-level-loop-checks/`.
No selected-mechanism build, full matrix, installed-binary audit, root/intermediate
selected-level application trajectory, recursive MPI EB integration, temporal
convergence or physical validation was run. Halo preparation still spans all
levels, requiring valid donor fields throughout the tree.
