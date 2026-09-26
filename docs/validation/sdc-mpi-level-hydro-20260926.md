# MPI SDC level hydro and temporal donors — 2026-09-26

Implementation `dabbb771`, base `691c9e36`; PeleF 0.246.0.
PeleC pin `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Implementation

The MPI SDC hydro stage accepts optional `target_depth`, packed parent state
start/end trees and fractions, and separate signed predictor-source start/end
trees and fractions. All ranks agree on selection and endpoint presence before
communication. Existing distributed ghost sampling validates interpolation
fractions collectively and fetches only requested owner-local endpoint cells.
Start endpoints replace ancestor fields in private copies. State endpoints have
`nvar+1` components including temperature; source endpoints have `nvar`.
Only ancestor endpoint payloads need allocation on their owners.

Inactive owned rates are zero; inactive and nonowned flux payloads are absent.
Inputs remain unchanged, and rates/fluxes publish only after collective success.
The existing common-step path remains the default. This stage does not itself
perform recursive stepping or coarse/fine synchronization.

## Focused evidence

Artifacts: `build/sdc-mpi-level-hydro-checks/` (local).
GNU Fortran 13.3.0, Open MPI 4.1.6, existing Ninja Debug tests-enabled build
`build/spray-amr-reactive-debug`; compiler/backend unchanged.

- Provenance passed: `python3 tools/check_mechanism_provenance.py --bundle
  mechanisms/h2o2_full.json`; SHA
  `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Build passed: `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
  build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_hydro_3d
  --parallel 2` (`build.log`). Initial build failed because the added test used
  an undeclared component-count name; replaced with `size(q)` and rebuilt.
  Initial failure retained in `initial-build-failure.log`. Warnings remain.
- CTest **2/2 passed**, 6.35 s total: `unit_mpi_amr_patch_tree_sdc_hydro_3d_4`
  1.95 s, `unit_mpi_amr_patch_tree_ppm_3d_2` 4.35 s. Command:
  `env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86
  OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir
  build/spray-amr-reactive-debug -R
  '^unit_mpi_amr_patch_tree_(sdc_hydro_3d_4|ppm_3d_2)$' --output-on-failure`.
- Extended the SDC test's nonuniform three-level tree, owners `[0,1,0]`,
  four ranks with two idle. Deepest-level hydro uses state fractions `[.25,.5]`
  and source fractions `[.75,.5]`, changed parent density/energy and signed
  source endpoints. Serial/MPI rates and all three integrated flux arrays match
  exactly. Inactive rates are zero and only the selected owner's flux exists.
- A differing state fraction on the last idle rank is rejected collectively;
  prior rates and x flux remain exact. Existing full-tree periodic/outflow
  parity and one-owner invalid-energy rollback checks also passed.
- Python syntax for `tools/check_project_contract.py`, project contract,
  touched Fortran line lengths and `git diff --check` passed.

No full matrix, selected build, application trajectory, recursive MPI SDC,
physical validation or installed audit was run. Temporal transport and the
coupled MPI level/recursive drivers remain to be connected.
