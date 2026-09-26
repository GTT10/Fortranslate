# MPI SDC level transport and temporal donors — 2026-09-26

Implementation `f674059c`, base `0ac142d3`; PeleF 0.246.0.
PeleC pin `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Implementation

The MPI SDC transport stage accepts optional `target_depth`, parent packed
state start/end trees and interpolation fractions. Packed states contain gas
components and temperature. Selection and endpoint presence are agreed across
ranks. Existing distributed ghost sampling validates fractions collectively
and fetches requested endpoint cells. Only owner-local ancestor endpoint
payloads need allocation. Start endpoints replace private ancestor copies.

Only the selected level contributes rates, integrated fluxes and the global
minimum limiter coefficient. Inactive owned rates are zero; inactive and
nonowned fluxes are absent. Inputs remain unchanged and outputs publish after
collective success. Common-step behavior remains the default.

## Focused evidence

Artifacts: `build/sdc-mpi-level-transport-checks/` (local). GNU Fortran 13.3.0,
Open MPI 4.1.6, existing Ninja Debug tests-enabled build
`build/spray-amr-reactive-debug`; compiler/backend unchanged.

- Provenance passed: `python3 tools/check_mechanism_provenance.py --bundle
  mechanisms/h2o2_full.json`; SHA
  `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Build passed: `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
  build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_transport_3d
  --parallel 2` (`build.log`). Compiler warnings remain.
- CTest **2/2 passed**, 10.97 s total: `unit_mpi_amr_patch_tree_sdc_transport_3d_4`
  1.74 s; `unit_mpi_amr_patch_tree_transport_3d_2` 9.18 s. Command:
  `env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86
  OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir
  build/spray-amr-reactive-debug -R
  '^unit_mpi_amr_patch_tree_(sdc_transport_3d_4|transport_3d_2)$' --output-on-failure`.
- Extended the nonuniform three-level SDC fixture, four ranks, owners `[0,1,0]`
  with two idle ranks. Deepest-level transport uses parent fractions `[.25,.5]`
  and end temperatures raised by 2 percent. Serial/MPI rates, all three flux
  arrays and limiter match exactly. Inactive rates and flux ownership are checked.
- A differing fraction on the last idle rank is collectively rejected; prior
  rates, x flux and limiter remain exact. Existing periodic and moving
  isothermal-wall stage parity and bad-coefficient rollback tests also passed.
- Python syntax for `tools/check_project_contract.py`, project contract,
  touched Fortran line lengths and `git diff --check` passed.

No application trajectory, selected build, recursive MPI SDC, full matrix,
installed audit or physical validation was performed. The coupled MPI level
and recursive drivers still need to use these stage interfaces.
