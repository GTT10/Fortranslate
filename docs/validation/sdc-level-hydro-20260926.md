# Level-selective SDC hydro and temporal donors — 2026-09-26

Implementation `6324431b`, base `63b7b39f`; PeleF 0.246.0.
Pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Implemented contract

`compute_sdc_amr_patch_tree_hydro_3d` accepts optional `target_depth`.
Unselected rate fields are zero and their integrated-flux components remain
unallocated. Default behavior still computes every patch.

Optional `coarse_start`, `coarse_end`, `coarse_fractions` provide packed conserved
state plus temperature endpoints for ancestors. Independently,
`source_start`, `source_end`, `source_fractions` supply signed predictor-source
endpoints. Each group must be complete and requires an explicit target depth.
Fractions have one entry per ancestor level and must be finite in [0,1].
Only ancestor fields need allocation in endpoint trees. The existing tree ghost
sampler handles interpolation; state composition limiting is retained, while
source fields use signed scalar interpolation. Outputs publish only on success.

This is a building block for recursive SDC, not a recursive time integrator.
The caller still owns parent intervals, level iteration barriers, source
estimates, transport/reaction stages, and reflux after child intervals.

## Completed focused checks

Artifacts: `build/sdc-level-hydro-checks/` (local).
GNU Fortran 13.3.0, Ninja Debug, existing tests-enabled
`build/spray-amr-reactive-debug`; compiler/backend unchanged.

- `python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json`
  passed, SHA `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug
  --target test_amr_patch_tree_hydro_3d --parallel 2` passed. Logs `build.log`
  and `build-test.log`. Existing compiler warnings remain.
- `unit_amr_patch_tree_ppm_3d`: **1/1 passed**, 11.89 s (11.94 s total).
  The periodic SDC checks now compare each selected level exactly with the
  full-tree stage, check inactive zero rates/absent fluxes, reproduce old
  temporal endpoints exactly, observe changed fine rates at half parent time,
  and verify rate rollback for an out-of-range fraction. Existing hydro checks
  remain. See `ctest.log`.
- Command: `env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH
  HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest
  --test-dir build/spray-amr-reactive-debug -R '^unit_amr_patch_tree_ppm_3d$'
  --output-on-failure`.
- Python syntax for `tools/check_project_contract.py`, project contract,
  Fortran line lengths and `git diff --check` passed.

No application run, selected rebuild, MPI check, recursive SDC solve,
installation audit, independent full-solver parity or physical validation was
performed. The temporal tests exercise periodic coarse/fine ghost sampling;
they do not establish temporal parity for all physical boundary conditions.
