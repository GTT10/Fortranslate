# Regular gas SDC reactor selection — 2026-09-27

Implementation `558ab778`, based on `b7b3e51b`, local main, version 0.246.0.

## Change

`pelef_reactive_3d` and `pelef_reactive_3d_selected` accept
`sdc_reactor='native'|'cvode'` in `reactive_3d`. Native remains the default.
CVODE requires reacting SDC and a SUNDIALS-enabled build. The frontend supplies
the forced CVODE callback explicitly; the shared driver verifies that input and
backend agree and prints the actual SDC reactor. Fixed and selected executables
link separate backend libraries against their corresponding core runtimes.

The configuration reader has an explicit CVODE capability flag, default false.
Other applications using that reader therefore reject the new selection until
their hierarchy paths are connected. Unknown values and nonreacting CVODE
requests also fail before output. Bundle `explicit`/`implicit` metadata is
unchanged; the new field selects the forced SDC reactor independently.

## Independently completed checks

GNU Fortran 13.3.0, Ninja Debug. SUNDIALS 7.2.0-enabled build
`build/forced-cvode-debug` retains the pinned H2O2 bundle/source configuration
from [the backend record](forced-cvode-20260927.md). The rejection check uses
SUNDIALS-disabled `build/spray-amr-reactive-debug`. Compiler/backend settings
were not changed within either build directory.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/forced-cvode-debug --target pelef_reactive_3d pelef_reactive_3d_selected --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/forced-cvode-debug
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/forced-cvode-debug -R '^application_sdc_cvode$' --output-on-failure
python3 -m py_compile tests/application/check_sdc_cvode_app.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Application CTest **1/1 passed**, 63.29 s, total 63.31 s, no failures/skips.
It runs fixed CVODE, selected CVODE, native default, unknown-selector rejection
and nonreacting-CVODE rejection. The matching-composition H2/O2 fixture is a
uniform periodic 4³ grid, 1100 K, 300000 Pa, molecular transport, nonzero
external force, two SDC iterations and final time 1e-8 s.
Fixed/selected CVODE CSV bytes agree exactly. Maximum native/CVODE temperature
difference is 4.3807667680084705e-6 K (criterion 1e-4 K); chemistry raises the
temperature by more than 0.001 K. The normalized Euler inventory residual is
1.6612134626384133e-17 (criterion 3e-12).

A separate subprocess runs the same fixed CVODE input with the disabled-build
executable, expecting `SDC CVODE backend is unavailable` and no output CSV.
That check passed; its input and log are retained in
`build/sdc-cvode-app-checks/unavailable/`. The checker also supports
`--unavailable PATH` to repeat that negative case alongside the others.

Builds, provenance, Python syntax, project contract and whitespace checks passed.
The initial compile exposed a duplicate configuration member/missing local
declaration; correcting declaration placement resolved it. Existing compiler
warnings remain. Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Logs are in `build/sdc-cvode-app-checks/`, including preserved `LastTest.log`.
Final per-case inputs/CSVs/logs:
`build/forced-cvode-debug/sdc-cvode-app/sdc_cvode_b7594iuk/`.

## Limits

AMR/EB/MPI/spray SDC reactor selection remains unconnected. These regular gas
applications have no restart interface. Nonuniform chemistry, other mechanisms,
CVODE with physical boundaries, installed binaries, full configuration matrix,
external full-solver parity and physical validation were not checked here.
