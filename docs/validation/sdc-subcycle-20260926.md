# Serial recursive AMR SDC — 2026-09-26

Implementation `3cfdc51d`, base `888bc720`; PeleF 0.246.0.
Pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Implemented behavior

`advance_amr_patch_tree_sdc_subcycle_3d` recursively advances all patches of
one level, then advances the next level with refinement-ratio substeps.
Parent state and predictor-source endpoints remain available throughout child
intervals. Integrated child fluxes accumulate before the existing depth-bounded
reflux/restriction operation. Temperatures recover after each synchronization.
External-source inventories accumulate only over selected-level leaf ownership.

Reaction estimates persist across substeps and are returned after the entire
coarse interval succeeds. Candidate state, temperatures, estimates, inventories,
reflux/transport diagnostics and optional int64 level counts remain private on
failure. Expected level counts are checked, including multiplication overflow.

The serial application now dispatches SDC with `subcycling=.true.`. Existing
subcycle timestep estimates and checkpoint context remain in use; SDC V3/V4
records carry the persistent source state. MPI and EB SDC are not enabled.

## Completed focused checks

Artifacts: `build/sdc-subcycle-checks/` (local).
GNU Fortran 13.3.0, existing Ninja Debug tests-enabled
`build/spray-amr-reactive-debug`; no compiler/backend change.

- Provenance passed: `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`;
  SHA `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Builds passed using `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
  build/spray-amr-reactive-debug --target TARGET --parallel 2`, for targets
  `test_amr_patch_tree_split_3d` and `pelef_amr_tree_reactive_3d`.
  Logs `build.log`, `build-test.log`, `build-app.log`; warnings remain.
- Focused CTest `unit_amr_patch_tree_sdc_3d`: **1/1 passed**, 17.36 s
  (17.40 s total), `ctest.log`. The added nonuniform three-level reacting case
  executes **1, 2, 4** intervals, includes molecular transport and forcing,
  closes Euler/element inventories within the existing 2e-10 relative limits,
  and has exactly restricted parent states after synchronization. An invalid
  hydro coefficient preserves state, temperature, diagnostic and count sentinels.
- CTest command uses unset DISPLAY/LD_LIBRARY_PATH/LIBRARY_PATH and
  `HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none`, with
  `ctest --test-dir build/spray-amr-reactive-debug -R '^unit_amr_patch_tree_sdc_3d$'
  --output-on-failure`.
- `python3 build/sdc-subcycle-checks/run_checks.py`: **three application runs
  passed** (continuous, checkpoint stop, restart). Fixed full H2/O2, two-level
  periodic grid, two SDC iterations, chemistry/transport/forcing, 3e-9 s final
  time. Continuous run uses two coarse steps. Indexed restart state and budget
  CSVs are byte-identical to continuous results (`comparison.log`). Inputs,
  runner, checkpoint, CSVs and logs retained.
- Maximum relative Euler inventory residual after external/boundary accounting:
  **7.535819449199924e-16** (`inventory.txt`). Minimum transport theta: 1.
- Python syntax for local runner and `tools/check_project_contract.py`,
  project contract, Fortran line lengths and `git diff --check` passed.

No selected rebuild/run, recursive regrid application, nonperiodic recursive
SDC case, chemistry-disabled recursive case, MPI/EB SDC, install audit,
full-solver pinned parity or physical validation was performed. The shared
selected driver is connected but not independently qualified in this change.
