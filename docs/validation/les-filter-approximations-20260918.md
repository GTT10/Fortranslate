# LES filter approximations — 2026-09-18

Implementation `240fd7763d9633b96d1170e033868f86723afbb7`, based on `7e006bc`.
The shared weight builder now adds `box3`, `box5`, `gaussian3`, `gaussian5`,
`box3opt`, `box5opt`, `gaussian3opt`, `gaussian5opt`. These implement the pinned
PelePhysics three/five-point polynomial and optimized approximation formulas.
Gaussian three-point polynomial weights equal box three-point weights. Optimized
ratio tables cover 1–10. Above ten, box optimized variants fall back to standard
box, Gaussian optimized variants to the corresponding polynomial approximation.
The standard box fallback rejects odd ratios, matching its reference constraint.

Reference: PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`,
`Source/Utility/Filter/Filter.cpp`; pinned SHA256 is in `references/les_filter.json`.
Negative weights are preserved. All serial filter consumers accept these names;
MPI dynamic coefficient and coupled-spray control IDs now use a shared mapping.
The existing API ratio cap is 64. Finite-value, normalization and EOS guards are
unchanged: wide signed approximations can be rejected. For example, direct
five-point formulas at ratio 64 have weight-sum roundoff exceeding the existing
`1e-12` absolute normalization gate. This increment does not loosen that gate or
claim positive-state preservation for signed filters.

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; unchanged configurations.
Fixed tests-enabled `build/spray-amr-reactive-debug` built `test_spray_core` and
`test_les_mpi_dynamic`. Selected H2O2 tests-disabled
`build/spray-amr-reactive-selected-debug` built `pelef_selected_spray_runtime`.

Focused CTest **4/4 passed**, 1.15 seconds, zero failures/skips. Existing standard
filter checks remain. New checks cover symmetry/normalization for all eight
approximations at ratios 1–10, signed periodic impulse conservation, polynomial
second/fourth moments, selected optimized table ratios, retained negative weights,
and pinned fallback behavior. MPI two/four-rank tests add box5opt/gaussian3opt
configurations to exact serial coefficient/RHS/diffusion/state/temperature parity,
conservation and collective rollback checks. Provenance, Python syntax for the
contract checker, project contract and whitespace checks passed. Selected runtime
was compiled, not exercised. No full matrix, external solver parity or physical
validation ran.

## Commands and artifacts

Logs: `build/les-filter-approximations-checks/{provenance,build,build-selected,ctest}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core test_les_mpi_dynamic --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_(filter|filter_approximations|mpi_dynamic_[24])$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
