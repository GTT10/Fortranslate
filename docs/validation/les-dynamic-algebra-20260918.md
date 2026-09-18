# Dynamic Smagorinsky algebra — 2026-09-18

## Implementation and reference

Implementation `8568f5697e11f0ed46e1b65634e0d44b5372f6bd`, based on `2ad098e`.
`les_dynamic_smagorinsky_mod` implements cell-local grid/test quantities and
Germano coefficient contractions from supplied derivatives and filtered inputs.
It returns density-weighted velocity products, velocity-temperature products,
deviatoric stress, isotropic strain quantity and thermal flux. Coefficients are
ordered `Cs^2`, `CI`, `Cs^2/PrT`; off-diagonal terms are contracted twice through
full symmetric matrices. CI is bounded above by one. Spatial differentiation,
spatial filtering, coefficient smoothing and dynamic flux/application dispatch
are not implemented by this increment.

Reference: pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`,
[LES.H](https://raw.githubusercontent.com/Pele-Suite/PeleC/bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8/Source/LES.H)
(`get_sfs_stresses_cc`, `pc_dynamic_smagorinsky_quantities`,
`pc_dynamic_smagorinsky_coeffs`) and
[Constants.H](https://raw.githubusercontent.com/Pele-Suite/PeleC/bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8/Source/Constants.H).
Exact source SHA256 values are in `references/les_dynamic.json` and downloaded
sources in `build/les-reference/`. The pinned `1e-8` denominator regularizer and
coefficient floor are retained; unit-rescaled dynamic LES equivalence is not
established. Invalid inputs return zero outputs and false status. Exactly zero
regularized CI denominator rejects; clipping before division prevents overflow
for very small nonzero denominators without changing the finite formula.

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug. Fixed tests-enabled build
`build/spray-amr-reactive-debug` and selected H2O2 tests-disabled build
`build/spray-amr-reactive-selected-debug` retain their existing configurations.

- Mechanism provenance passed.
- Fixed `test_spray_core` and selected `pelef_selected_spray_runtime` built.
- Final focused CTest **2/2 passed**, 0.04 seconds, zero failures/skips:
  existing Smagorinsky algebra plus dynamic quantities/contractions.
- Dynamic checks cover analytic shear, filter-width scaling, weighted products,
  tensor symmetry, nonzero off-diagonal contractions, signed heat contraction,
  unclipped and upper-clipped CI, zero/negative coefficient floors, singular CI,
  nonfinite input and invalid filter ratio rejection.
- Python syntax for the contract checker, project contract and diff whitespace
  checks passed. No full matrix, MPI dynamic runtime, installed-binary audit,
  external solver parity or physical validation was run.

The first focused run passed 2/2 in 0.05 seconds; the test was then extended with
an interior CI value and rerun as above. No failed builds or test runs occurred.

## Commands and artifacts

Logs are under `build/les-dynamic-algebra-checks/`: `provenance.log`, `build.log`,
`build-final.log`, `build-selected.log`, `build-test-final.log`, `ctest.log` and
`ctest-final.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_(les_algebra|les_dynamic_algebra)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
