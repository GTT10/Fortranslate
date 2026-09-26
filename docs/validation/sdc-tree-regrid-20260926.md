# Serial AMR SDC regridding — 2026-09-26

Implementation `4148ce64`, base `e57d4b58`; PeleF 0.246.0.
Pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Implementation

The reactive tree remap accepts paired optional old/new reaction-source trees
with `nspecies+1` components. Signed rates use the existing geometric overlap
algorithm, without positivity clipping. It remaps leaf values and averages
covered parents. Gas, recovered temperature and reaction estimates publish
only after all stages succeed. The driver forwards its persistent reaction
tree for scheduled and temperature-triggered regrids. The source inventory
ledger remains a cumulative domain inventory and does not need remapping.
Recursive SDC timestepping remains unsupported.

## Focused evidence

Artifacts: `build/sdc-tree-regrid-checks/` (local).
GNU Fortran 13.3.0, existing Ninja Debug tests-enabled build
`build/spray-amr-reactive-debug`; compiler/backend unchanged.

- Provenance passed: `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`;
  SHA `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Build passed: `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
  build/spray-amr-reactive-debug --target pelef_amr_tree_reactive_3d
  test_amr_patch_tree_remap_3d --parallel 2` (`build.log`). Warnings remain.
- `unit_amr_patch_tree_remap_3d`: **1/1 passed**, 0.01 s (0.05 s total),
  including signed nonuniform reaction-source conservation across changed
  branches/refinement ratios and exact output retention after an EOS failure.
  Existing dense overlap oracle and gas tests remain. See `ctest.log`.
- Application runner: `python3 build/sdc-tree-regrid-checks/run_checks.py`.
  **Four successful runs**: scheduled continuous, checkpoint stop, restart,
  automatic adaptation. Each uses fixed H2/O2, reactions, transport, forcing
  and two SDC iterations through 1e-9 s.
- Scheduled case removes refinement at 5e-10 s and reinstates it at 7.5e-10 s.
  Continuous run finishes in three steps. Restart from the checkpoint after
  the first remap produces byte-identical state and budget CSVs. Maximum gas
  remap inventory error reported as zero. `comparison.log`, all inputs,
  checkpoints, CSVs and logs retained.
- Automatic case removes refinement at step 1 (8.5781409157328871e-10 s),
  then completes at step 2. Two adaptation checks, one topology change.
- Initial application attempt failed before stepping because the runner wrote
  unquoted absolute paths in a Fortran list-directed schedule. Quoting those
  paths corrected the input. Initial failure retained as `initial-path-failure.log`;
  it is excluded from successful counts.
- CTest used `env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH
  HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest
  --test-dir build/spray-amr-reactive-debug -R '^unit_amr_patch_tree_remap_3d$'
  --output-on-failure`.
- Python syntax for the local runner and `tools/check_project_contract.py`,
  project contract, Fortran line lengths and `git diff --check` passed.

The application fixture is spatially uniform; the unit test supplies nonuniform
signed rate coverage. No selected rebuild/run, MPI test, recursive SDC case,
physical validation, installed audit or pinned full-solver parity was run.
The existing piecewise volume-overlap algorithm is not a claim of matching
all AMReX interpolation choices.
