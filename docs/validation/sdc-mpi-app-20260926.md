# MPI AMR SDC application — 2026-09-26

Implementation `8784880e`, base `4e8b3503`; PeleF 0.246.0.
PeleC pin `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Implementation

Fixed/selected MPI frontends accept common-step SDC and external forcing.
The driver retains owner-local reaction estimates and the cumulative external
inventory. Checkpoint save/restore forwards both fields and verifies their
presence against the configuration. Restart conservation includes external
forcing. Scheduled/automatic remapping transfers signed estimates with generic
volume overlap before publishing gas, temperature and estimates. Rebalancing
packs gas and estimates into one migration payload and publishes the new owner
map with all fields. MPI SDC subcycling is explicitly rejected at input.

## Focused evidence

Artifacts: `build/sdc-mpi-app-checks/` (local). GNU Fortran 13.3.0,
Ninja Debug, Open MPI 4.1.6. Existing fixed tests-enabled and selected
tests-disabled builds; compiler/backend unchanged.

- Provenance preflight passed with `python3 tools/check_mechanism_provenance.py
  --bundle mechanisms/h2o2_full.json`; SHA
  `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Builds passed using `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
  BUILD --target TARGET --parallel 2`: fixed BUILD
  `build/spray-amr-reactive-debug`, TARGET `pelef_mpi_amr_tree_reactive_3d`;
  selected BUILD `build/spray-amr-reactive-selected-debug`, TARGET
  `pelef_mpi_amr_tree_reactive_3d_selected`. Logs `build.log` and
  `selected-build.log`. Compiler warnings remain.
- `python3 build/sdc-mpi-app-checks/run_checks.py`: **3/3 application runs
  passed**, using `/usr/bin/mpiexec`: four-rank continuous, four-rank checkpoint
  stop, two-rank restart. Fixed H2/O2, two SDC iterations, chemistry, molecular
  transport and forcing; root/child grids each 4 cubed. Coarsening at 5e-10 s,
  refinement at 7.5e-10 s, final time 1e-9 s in three steps. The indexed
  checkpoint is saved after the first remap. Final state and budget CSVs match
  byte for byte between continuous and changed-rank restart (`comparison.log`).
- Maximum reported remap error is zero. First five conserved components have
  maximum absolute ledger residual divided by initial magnitude
  `1.5088880917347677e-16`. CSVs also match the existing serial artifacts in
  `build/sdc-tree-regrid-checks/`; that serial calculation was not rerun here.
- Python syntax for the runner and `tools/check_project_contract.py`, project
  contract and `git diff --check` passed. No separate unit suite rerun.

The application fixture is spatially uniform. No selected application run,
actual measured-cost migration, automatic adaptation, nonperiodic application,
recursive MPI SDC, physical validation, full matrix or installed audit was
performed. The source migration path is implemented but actual rebalancing
needs a dedicated case. Earlier stage/unit evidence remains separate.
