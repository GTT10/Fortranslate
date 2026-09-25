# Single-domain EB spray viscosity — 2026-09-25

Code `7a36de2`, parent `e799e29`. Fixed/selected EB spray input now enables
finite nonnegative `difmag`. Both allocated and fixed-particle coupling APIs
forward it to the gas full step. Validation occurs before private particle
halves; candidate state/particle publication remains transactional. Both
application wall-impact branches forward the coefficient.

Checkpoint context uses explicit legacy configuration components and appends
`difmag-v1` only for positive coefficients. This also repairs the unintended
zero-context format change caused by adding a field to the configuration type
in `eb5987b`. Interim checkpoints containing that extra zero field are not
claimed compatible. The older component sequence is retained.

## Completed checks

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug; reused
`build/spray-amr-reactive-debug`. Fixed executable built. Focused CTests **2/2
passed**, positive input/restart 1.40 s and existing force/restart 1.09 s,
2.53 s total. No failures or skips; compiler warnings remain.

Coefficient 0.037 with particle acceleration, chemistry and transport gives
exact uninterrupted/restarted gas, particle and film CSVs and printed ledgers.
Changed coefficient, negative and NaN inputs reject and preserve those outputs.
Existing context/geometry/corruption/publication checks also run. A separate
GNU Fortran serialization probe compares the configuration type from `d948f3c`
against the explicit legacy component sequence: default output bytes match.
This is a format probe, not a historical executable restart audit.

Selected frontend, positive wall-impact dynamics, MPI EB spray and the full
configuration matrix were not run. This does not establish nozzle, diesel
ignition or physical validation. MPI EB spray input remains disabled for
positive coefficients pending its coupling connection.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_(difmag|force)_application$' --output-on-failure
python3 -m py_compile tests/application/check_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

All listed checks passed. Artifacts in `build/eb-spray-difmag-checks/`:
`build.log`, `LastTest.log`, `context_format.f90`, `context_format.log`,
`legacy.txt`, `explicit.txt`. The probe was compiled with
`gfortran context_format.f90 -o context_format` and executed in that directory.
