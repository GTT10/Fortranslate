# Regular 3D pressure outlet application — 2026-09-18

Implementation `2bb6bba`, base `16b2cb0`.

The shared `reactive_3d` namelist now accepts `outlet_pressure(6)`, initialized
to zero. Positive finite entries up to 1e15 Pa require an outflow face.
The fixed and selected regular 3D front ends explicitly enable this input and
forward it through the full reactive wrapper. Other shared-reader applications
reject nonzero targets until their numerical paths are connected; they cannot
silently ignore the new input. The regular application has no checkpoint path.

The added unit checks cover input parsing, reader capability rejection,
positive-pressure face validation and negative-pressure rejection, alongside
the existing six-face ghost/hydro checks. The application fixture enables
chemistry and molecular transport with characteristic PLM/HLLC and a 102325 Pa
x-high target from initially uniform 101325 Pa gas.

GNU Fortran 13.3.0 / Open MPI 4.1.6, Debug/Ninja; existing fixed tests-enabled
and selected H2O2 tests-disabled configurations. Selected executable compiled;
no selected runtime or MPI runtime qualification in this increment.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_3d test_reactive_boundaries_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_3d_selected --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(unit_reactive_boundaries_3d|reactive_boundary_3d_pressure_outlet)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/pressure-outlet-app-checks/`. AMR/MPI pressure connections,
characteristic inflow, acoustic reflection and physical validation remain open.

## Results

CTest **2/2 passed**, zero skips, 98.67 s (unit 5.86 s, application 92.77 s).
Both builds, provenance, Python syntax, project contract and diff checks passed.
The unit target was rebuilt after adding the reader round-trip checks, before
CTest. No failed or interrupted runs; no tolerance changes.

The application completed one step to 2e-7 s. Maximum speed was
0.3169704607 m/s, x-momentum change -2.0206212598e-8 kg m/s, and the maximum
boundary-accounted normalized Euler conservation residual was 3.7337212315e-18
(existing gate 5e-11, normalization floor 1). Species integral change was
2.0411317908e-11 kg; this includes chemistry and boundary exchange.
These are execution/conservation observations, not acoustic nonreflection or
physical validation.

Logs: `provenance.log`, `build.log`, `build-final.log`, `build-selected.log`,
`ctest.log` in the artifact directory. The detailed CTest log is copied as
`LastTest.log`; output CSV is
`build/spray-amr-reactive-debug/tests/boundary_3d_pressure_outlet/pressure_outlet.csv`.
