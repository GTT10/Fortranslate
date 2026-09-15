# Local spray integration, 2026-09-15

Pulled the checkpoint branch, updated local main to integrated main
`5a2c818899f7213dad2fbb3c6808558711aced10`, then fast-forwarded the existing
spray implementation to `100e21e`. No remote branch was changed.

The first clean GNU build failed: substituting the 64-character mechanism
hash into `app/pelef_spray_3d.F90.in` exceeded Fortran's 132-column limit.
Both fixed and selected application templates now continue the hash argument
on a separate line. The tested source is `100e21e` plus those two template
changes; this record adds no numerical change.

## Focused verification

GNU Fortran 13.3.0 (`/usr/bin/gfortran`), Ninja, Python 3.12.3;
the optional configuration uses existing SUNDIALS 7.2.0 official Fortran
libraries and Cantera 3.2.0. Both build directories were newly configured.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake -S . -B build/port-spray-debug -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran \
  -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_TESTS=ON
cmake --build build/port-spray-debug --parallel 8 \
  --target pelef_spray_3d test_spray_core test_spray_integration
ctest --test-dir build/port-spray-debug -L spray --output-on-failure --no-tests=error

build/cantera-venv/bin/python tools/prepare_spray_reference.py \
  --output-dir build/port-methanol-reference
cmake -S . -B build/port-spray-cvode -G Ninja \
  -DCMAKE_BUILD_TYPE=Release -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran \
  -DPELEF_ENABLE_SPRAY=ON -DPELEF_ENABLE_TESTS=ON \
  -DPELEF_ENABLE_SUNDIALS=ON -DPELEF_ENABLE_CANTERA_REFERENCE=ON \
  -DPython3_EXECUTABLE="$PWD/build/cantera-venv/bin/python" \
  -DCMAKE_PREFIX_PATH="$PWD/build/sundials-7.2-fortran/install" \
  -DPELEF_MECHANISM_BUNDLE="$PWD/build/port-methanol-reference/ffcm1_red.json" \
  -DPELEF_MECHANISM_SOURCE="$PWD/build/port-methanol-reference/ffcm1_red.yaml"
cmake --build build/port-spray-cvode --parallel 8 \
  --target pelef_spray_3d_selected pelef_spray_3d test_spray_integration
ctest --test-dir build/port-spray-cvode --parallel 3 --output-on-failure \
  --no-tests=error \
  -R '^(spray_integration_cvode|spray_cantera_methanol-chemistry|spray_cantera_methanol-spray|spray_reactive_cone_restart)$'
python3 tools/check_project_contract.py --source .
git diff --check
```

- Debug spray suite: **24/24 passed**, 61.87 seconds.
- Release CVODE/methanol/cone selection: **4/4 passed**, 28.08 seconds.
- Methanol gas chemistry: maximum Cantera history temperature difference
  **0.0301034 K**; maximum final mass-fraction difference **1.049e-7**.
- Closed methanol spray: all liquid evaporated by 50 microseconds;
  relative system mass/energy residuals **1.525e-15 / 8.140e-16**.
- Provenance, project contract, diff whitespace and syntax of the three
  preparation/chemistry/cone Python tools passed.

Local artifacts: `build/port-spray-{debug,cvode}/configure.log`, `build.log`,
`focused-tests.log`, and `Testing/Temporary/LastTest.log` (individual outputs).
The pinned methanol source, license and generated bundle are under
`build/port-methanol-reference/`.

This is focused integration evidence, not a full suite or physical spray
validation. General 3D boundaries/AMR/EB, distributed particles and the other
open capabilities in the completion roadmap remain unfinished.
