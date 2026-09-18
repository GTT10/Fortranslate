# Dynamic LES spray coupling — 2026-09-18

Implementation `fce166db7ba4387abd5101d6e231d7870db2d74b`, based on `c2406e0`.
The serial `advance_spray_reactive_3d` library API now accepts
`les_model='dynamic'` and trailing optional `les_test_filter`, `les_test_ratio`,
`les_grid_ratio`. Defaults are box, two, one. Both SGS half steps use the dynamic
SSPRK2 operator within the existing candidate gas/particle/film transaction.
Dynamic LES requires periodic faces, including when disabled. Filter validity
is checked before stepping. Boundary exchange from periodic dynamic SGS is zero.

`cs=0` keeps the existing disabled convention. For dynamic LES only, positive
`cs` enables the model but its magnitude does not set coefficients; fixed `prt`
and `sct` likewise do not set dynamically calculated heat transport. Dynamic SGS
has no species flux. Existing numeric ranges remain validated. Default, WALE and
Vreman dispatch retains its previous stepping path. Application input/restart
controls and MPI dynamic coupling are not connected by this increment.
The frozen-coefficient timestep estimate retains the limitations documented in
[the stepping record](les-dynamic-step-20260918.md).

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug. Unchanged fixed tests-enabled
`build/spray-amr-reactive-debug` built `test_spray_integration`; selected H2O2
tests-disabled `build/spray-amr-reactive-selected-debug` built
`pelef_selected_spray_runtime`.

Focused CTest **2/2 passed**, 9.80 seconds, zero failures/skips: existing
`spray_integration_coupled` and new `spray_integration_coupled_dynamic`.
The 5-cubed periodic dynamic fixture uses crossed gas velocities, implicit
chemistry, molecular transport, drag/heating/evaporation and a `2e-7` step.
Existing mass/energy relative `1e-12` and momentum absolute `1e-18` gates pass,
with active evaporation and chemistry. Nonperiodic configuration and invalid
filter ratio reject before mutation even at zero cs. A Gaussian dynamic half
step followed by invalid chemistry preserves gas/temperature/particles.
The existing invalid-model rollback check remains.

Provenance, Python syntax for the contract checker, project contract and
whitespace checks passed. No build/test failures occurred. The first CTest
command used a wrong name pattern and found zero tests; it is not counted as
passing evidence. The corrected exact pattern selected the two tests above.
Selected runtime was compiled, not exercised. No full matrix, MPI dynamic
runtime, external solver parity or physical validation was run.

## Commands and artifacts

Logs: `build/les-dynamic-coupled-checks/{provenance,build,build-selected,ctest,ctest-final}.log`.
The initial `ctest.log` contains the zero-test attempt; `ctest-final.log` contains
the completed two-test result.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_integration --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_integration_(coupled|coupled_dynamic)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
