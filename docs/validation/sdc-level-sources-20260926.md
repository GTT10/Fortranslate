# Level-local SDC transport and chemistry — 2026-09-26

Implementation `eab4093f`, base `6547f422`; PeleF 0.246.0.
Pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Contract

`compute_sdc_amr_patch_tree_transport_3d` accepts optional `target_depth` and
packed parent `coarse_start`, `coarse_end`, `coarse_fractions`. Only the selected
level produces rates and integrated flux; inactive rates are zero and flux
components are unallocated. The endpoint group must be complete, with finite
fractions in [0,1], one per ancestor level. Only ancestor endpoint fields need
allocation. The existing viscous tree ghost sampler performs interpolation.

`react_sdc_amr_patch_tree_3d` accepts optional `target_depth`. With selection,
output state/temperature/reaction trees must already be valid. Inactive fields
retain their output values, even if those differ from the interval's old state.
Selected patches integrate from the supplied old state and nonreacting rates.
This accommodates levels at different times. Without selection the existing
whole-tree allocation/behavior is preserved. Both APIs publish only on success.

These are stage APIs. Recursive interval scheduling, coupled level SDC
iterations, and synchronization after child steps are still required.

## Completed focused checks

Artifacts: `build/sdc-level-sources-checks/` (local).
GNU Fortran 13.3.0, Ninja Debug, existing tests-enabled
`build/spray-amr-reactive-debug`; no compiler/backend change.

- Provenance passed: `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`;
  SHA `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Build passed: `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
  build/spray-amr-reactive-debug --target test_amr_patch_tree_transport_3d
  test_amr_patch_tree_chemistry_3d --parallel 2` (`build.log`). Warnings remain.
- Focused CTest: **2/2 passed**, 66.77 s total with `-j2`.
  `unit_amr_patch_tree_transport_3d`: 34.58 s;
  `unit_amr_patch_tree_chemistry_3d`: 66.73 s. See `ctest.log`.
- Transport checks compare selected rates/fluxes exactly with the whole-tree
  result, assert inactive zero rates/absent fluxes, reproduce the old coarse
  endpoint exactly, and observe changed fine transport at half parent time.
- Reaction checks reproduce selected state, temperature and reaction rates
  exactly while retaining deliberately different inactive fields. Existing
  independent uniform-cell reaction comparisons and late-patch rollback remain.
- Command: `env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH
  HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest
  --test-dir build/spray-amr-reactive-debug
  -R '^unit_amr_patch_tree_(transport|chemistry)_3d$' --output-on-failure -j2`.
- Python syntax for `tools/check_project_contract.py`, project contract,
  Fortran line lengths and `git diff --check` passed.

No application run, selected rebuild, MPI test, recursive SDC solve, install
audit, full-solver reference parity or physical validation was performed.
Temporal transport additions were exercised in periodic cases.
