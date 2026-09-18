# MPI spray pressure outlet — 2026-09-18

Implementation `8270cee2f24147c281dadda97af0c3bda568d6a0`, based on `1c08eed`.

MPI hydro and molecular-transport wrappers, mixed reactive-spray coupling and
fixed/selected MPI application namelists accept trailing `outlet_pressure(6)`.
The default is zero. Finite values in `[0,1e15]` Pa and outflow-face assignment
are validated collectively. Exact pressure agreement is checked before updates;
nonzero pressure requires physical faces and, in coupling, the hydro path.
All gas operators receive the same pressure array. Nonzero pressure enters MPI
checkpoint identity; changing/removing it rejects restart. Parcel boundary
semantics and dynamic LES's periodic-only restriction remain unchanged.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing configurations unchanged.
Fixed tests-enabled `pelef_mpi_spray_3d` and selected H2O2 tests-disabled
`pelef_mpi_spray_3d_selected` built. Selected executable was not run.
Two distinct focused tests passed:

- `spray_mpi_isotropic_restart`: 6.98 s on the initial run.
- `spray_mpi_pressure_outlet_restart`: 11.18 s on the corrected-fixture rerun
  (CTest total 11.22 s, 1/1 passed).

The pressure fixture uses all-outflow faces and 320000 Pa on x-high, with
reacting mixed spray, molecular transport, WALE coefficient 0.17 and CI 0.1.
Checks cover exact same-rank restart, 2-to-4-rank continuation with the existing
`2e-11` relative/`1e-25` absolute CSV tolerance, mass/energy budgets, changed/
removed pressure identity rejection, negative pressure input rejection and
existing CI/model guards. No tolerances were changed.

Initial combined CTest was 1/2 passed: the new fixture used 100000 Pa against
its roughly 323000 Pa initial gas pressure, making the characteristic ghost
extrapolation inadmissible. Coupling correctly rejected all retries. The fixture
outlet was corrected to 320000 Pa; only the failed test was rerun. Initial logs
are retained. This does not qualify large pressure jumps or physical spray flow.
Provenance, Python syntax, project contract and whitespace checks passed.
No full matrix, external solver parity or experimental validation ran.

## Commands and artifacts

Logs: `build/spray-mpi-pressure-outlet-checks/`
(`provenance.log`, `build.log`, `build-selected.log`, `ctest.log`,
`ctest-low-pressure.log`, `ctest-pressure.log`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(pressure_outlet_restart|isotropic_restart)$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_pressure_outlet_restart$' --output-on-failure
python3 -m py_compile tests/spray/test_mpi_spray_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
