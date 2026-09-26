# Coupled level SDC interval — 2026-09-26

Implementation `dc266dda`, base `878b25a9`; PeleF 0.246.0.
Pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Contract

`advance_amr_patch_tree_sdc_3d` now optionally selects `target_depth`.
Level mode requires `integrated_flux` output and leaves restriction/reflux to
its caller. Other levels retain their state and temperature. Only the selected
level's leaf contribution enters the external-source inventory; physical root
flux inventory is returned at depth zero, with later fine correction left to
synchronization. Whole-tree mode retains its existing synchronization.

Fine-level calls require packed parent state endpoints, interval start/end
fractions, and parent predictor-source endpoints. Hydro samples parent state
at interval start and predictor sources at the midpoint. Old/new transport
samples parent state at interval start/end respectively. The existing level
reaction stage completes each SDC iteration without altering inactive levels.

Optional `old_predictor`/`new_predictor` outputs provide external plus transport
plus reaction source estimates for child prediction. Optional flux export is
hydro plus half old/new transport, integrated over this interval. State, source
estimates, diagnostics and fluxes publish only after successful completion.
This API advances one interval; it does not yet schedule recursive substeps.

## Completed focused checks

Artifacts: `build/sdc-level-loop-checks/` (local).
GNU Fortran 13.3.0, Ninja Debug, existing tests-enabled
`build/spray-amr-reactive-debug`; compiler/backend unchanged.

- Provenance passed: `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`;
  SHA `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Build: `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug
  --target test_amr_patch_tree_split_3d --parallel 2`. Final log `build-final.log`.
  Compiler warnings remain.
- Final focused CTest `unit_amr_patch_tree_sdc_3d`: **1/1 passed**, 11.86 s
  (11.91 s total), `ctest-final.log`. Earlier run also passed (11.66 s), before
  adding finite-output guards on optional predictor exports; it is not an
  additional distinct test case.
- Added checks advance the root level with chemistry, transport and forcing,
  compare its uniform result to the regular SDC reference within 1e-11 scaled
  tolerance, and assert unchanged finer state/temperature and absent finer flux.
  A coupled fine half-interval then consumes the returned parent state/source
  endpoints and leaves all other levels exactly unchanged. Fine interval has
  no root boundary ledger or reflux correction. Existing nonuniform whole-tree
  inventory, element, restriction and failure tests remain.
- Command: `env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH
  HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest
  --test-dir build/spray-amr-reactive-debug -R '^unit_amr_patch_tree_sdc_3d$'
  --output-on-failure`.
- Python syntax for `tools/check_project_contract.py`, project contract,
  Fortran line lengths and `git diff --check` passed.

No full recursive integration/synchronization case, application run, selected
rebuild, MPI run, install audit, full-solver reference parity or physical
validation was performed. Parent temporal coupling was exercised periodically.
