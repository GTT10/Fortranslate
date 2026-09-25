# SDC hierarchy hydro stage — 2026-09-26

Implementation `f66fa1ee`, parent `239ff81c`, local main, version 0.246.0.

`compute_sdc_amr_patch_tree_hydro_3d` builds PPM/PeleC hydro rates across all
patches from a shared old-state tree and predictor-source tree. State halos
use existing conservative coarse interpolation with species treatment and
physical boundary handling. Signed source fields use scalar MC interpolation
without a composition constraint; physical source halos extrapolate. The
stage returns cell rates and dt-integrated full-face fluxes transactionally,
without modifying state or applying reflux/restriction prematurely. These
fluxes use the existing AMR register layout and still need face-area weights
in the conservation correction.

## Completed checks

- GNU Fortran 13.3.0 Debug fixed tests-enabled build
  `build/spray-amr-reactive-debug`, target `test_amr_patch_tree_hydro_3d`, passed.
- `unit_amr_patch_tree_ppm_3d`: **1/1 passed**, 10.90 s test / 10.95 s total.
  Added root-only and branched multilevel periodic stage checks. Reconstructing
  the update from returned rates and applying existing tree reflux matches
  the existing PPM/PeleC hierarchy step within 1e-12 scaled state error.
  A nonzero momentum predictor source changes the computed hydro rate,
  exercising signed source interpolation. Existing boundary/CFL checks remain.
- Provenance, Python syntax, project contract and whitespace passed. No
  failed/skipped tests; existing compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_amr_patch_tree_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_amr_patch_tree_ppm_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/sdc-tree-hydro-checks/` contains build, provenance, CTest and
LastTest logs. Full hierarchy source/transport iteration and reaction coupling
remain to be orchestrated; the test applies reflux only as an isolated hydro
oracle. No complete AMR SDC, MPI or external solver equivalence is claimed.
