# Serial AMR SDC checkpoint/restart — 2026-09-26

Implementation `63cd9751`, base `1fcf715d`; PeleF 0.246.0.
Reference PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Change

V3 plain and V4 indexed records add an external-source ledger and optional
signed reaction-source fields. Gas record offsets retain their existing
meaning; reaction fields follow all gas records. V1/V2 writing is retained for
snapshots without the SDC ledger. The serial reader validates the complete
candidate before publishing it. Legacy header callers (including MPI) reject
V3/V4 unless they explicitly opt in. No MPI SDC continuation is claimed.

The serial application saves/restores reaction estimates and source inventory,
checks their presence against the configured integrator/chemistry, and includes
external sources in restart conservation checks. Configuration identity already
captures SDC iterations and forcing through the complete canonical config.

## Completed checks

Artifacts: `build/sdc-tree-restart-checks/` (local).
GNU Fortran 13.3.0, Ninja Debug, existing tests-enabled
`build/spray-amr-reactive-debug`; no compiler/backend changes.

- Provenance: `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`
  passed, SHA `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Build: `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug
  --target pelef_amr_tree_reactive_3d test_amr_patch_tree_checkpoint_3d --parallel 2`.
  Logs `build.log`, `build-test.log`; compiler warnings remain.
- Focused CTest `unit_amr_patch_tree_checkpoint_3d`: **1/1 passed**, 0.10 s
  (0.15 s total). Uses the existing test plus signed reaction estimates,
  exact V3/V4 round trips, legacy header rejection, context rejection with
  output retention, and nonreacting SDC records. Existing invalid-write and
  truncated V1 checks remain. Log `ctest.log`.
- `python3 build/sdc-tree-restart-checks/run_checks.py`: **five application runs
  passed** (one continuous; stop and restart for plain/indexed). Two-level
  periodic fixed H2/O2 case, two SDC iterations, chemistry/transport/forcing,
  two steps to 1e-9 s. Both restarted state CSVs and budget CSVs are byte-identical
  to continuous output. Exact inputs, runner, checkpoints, CSVs and logs retained;
  summary `comparison.log`. Based on committed `sdc.nml`/`sdc.dat`.
- CTest command used unset DISPLAY/LD_LIBRARY_PATH/LIBRARY_PATH and
  `HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none`, with
  `ctest --test-dir build/spray-amr-reactive-debug -R '^unit_amr_patch_tree_checkpoint_3d$'
  --output-on-failure`.
- Python syntax (`tools/check_project_contract.py` and local runner),
  project contract and `git diff --check` passed.

No selected-mechanism rebuild/application run, MPI suite, install audit, dynamic
SDC regrid, recursive SDC stepping, full solver parity or physical validation
was performed. Selected serial code shares this driver, but that fact is not
an independently rerun selected restart result.
