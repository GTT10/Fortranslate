# Recursive MPI SDC and application restart — 2026-09-27

Implementation `5ff33363`, base `5cd0282f`; PeleF 0.246.0.
PeleC pin `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Implementation

`advance_mpi_amr_tree_sdc_subcycle_3d` advances owner-local patch trees with
refinement-ratio child intervals. Parent state and predictor endpoints feed the
coupled level SDC stage. Child integrated fluxes accumulate before distributed
reflux/restriction, and temperature recovery follows synchronization. External
inventory accumulates over active leaves and local intervals. Root boundary
flux and reflux corrections form the boundary ledger. Optional 64-bit level
counts publish only after all expected intervals complete. State, temperature,
reaction estimates, diagnostics and counts remain private until full success.

The MPI application dispatches SDC `subcycling=.true.` to this driver. Existing
checkpoint context binds the subcycling mode and stores persistent reaction
estimates and the forcing ledger. The default common-step path is retained.

## Focused evidence

Artifacts: `build/sdc-mpi-subcycle-checks/` (local). GNU Fortran 13.3.0,
Open MPI 4.1.6, existing Ninja Debug tests-enabled build
`build/spray-amr-reactive-debug`; compiler/backend unchanged.

- Provenance passed: `python3 tools/check_mechanism_provenance.py --bundle
  mechanisms/h2o2_full.json`; SHA
  `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Builds passed: `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
  build/spray-amr-reactive-debug --target TARGET --parallel 2`, TARGET
  `test_mpi_amr_patch_tree_split_3d` (`build.log`) and
  `pelef_mpi_amr_tree_reactive_3d` (`app-build.log`). Warnings remain.
- CTest **1/1 passed**, 20.52 s (20.57 s total):
  `env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86
  OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir
  build/spray-amr-reactive-debug -R '^unit_mpi_amr_patch_tree_sdc_3d_4$'
  --output-on-failure` (`ctest.log`).
- Nonuniform three-level fixture, refinement ratios 2 then 3, owners `[0,1,0]`,
  four ranks with two idle. Two SDC iterations, reactions, all transport flags,
  forcing and dt=1e-10 s. Recursive state, temperature, reaction estimates,
  inventories, diagnostics and level counts match serial exactly. Invalid
  artificial viscosity on an idle rank rejects the interval after transport
  evaluation and retains state, temperature, reaction estimates, ledgers and
  counts. Existing common-step and level-stage checks also passed.
- `python3 build/sdc-mpi-subcycle-checks/run_checks.py`: **3/3 application runs
  passed** using `/usr/bin/mpiexec`: four-rank continuous, four-rank checkpoint
  stop, two-rank restart. Fixed full H2/O2, two-level uniform fixture, chemistry,
  transport, forcing, two SDC iterations and recursive stepping. Final time
  3e-9 s in two coarse steps; indexed checkpoint at step 1, 1.7156281831465774e-9 s.
  Continuous/restart state and budget CSVs match byte for byte (`comparison.log`).
- First five components' maximum ledger residual divided by initial magnitude:
  `7.535819449199924e-16`. CSVs also match existing serial artifacts in
  `build/sdc-subcycle-checks/`; that application was not rerun here.
- Python syntax for the local runner and `tools/check_project_contract.py`,
  project contract, touched Fortran line lengths and `git diff --check` passed.

No selected build/run, recursive regrid/rebalance trajectory, nonperiodic
recursive case, full matrix, physical validation or installed audit was run.
The full port remains incomplete, including EB SDC.
