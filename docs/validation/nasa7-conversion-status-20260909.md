# NASA7 dimensional-conversion failure status — 2026-09-09 JST

## Defect and bounded repair

Review of the integration candidate `f7f87b9af268f75995392bdaeb77206509123cce`
found that `nasa7_molar_properties` passed its `ok` output directly to
`nasa7_dimensionless_properties`. After successful dimensionless evaluation,
some dimensional-conversion overflow guards returned without resetting `ok`.
The routine therefore reported success with zero initialized outputs. The
mass-specific wrapper could propagate the same false success.

The repair resets `ok = .false.` immediately after successful dimensionless
evaluation, before any dimensional-conversion guard. The existing final
`ok = .true.` is retained after all validated outputs have been assigned.
The formulas, physical coefficients, polynomial limits and numerical
acceptance tolerances are unchanged. This is an error-propagation correction,
not a new thermodynamic model. Unlike the preceding CMake layout/CI commits,
it intentionally changes numerical-library failure behavior.

## Regression and red/green evidence

`test_nasa7_thermo` now constructs three finite synthetic coefficient sets at
1 K with valid synthetic temperature bounds. The dimensionless evaluation
must first succeed. Coefficients 1, 6 and 7 then independently force the
heat-capacity, enthalpy and entropy dimensional-conversion guards. Both the
molar and mass-specific interfaces must return failure and zero outputs.
The fixtures exercise rejection paths; they are not physical H2 data.

The new test with the original implementation failed with exit code 8 from
CTest and the diagnostic:

```
NASA7 molar conversion reported success after overflow rejection
```

After the status correction, the same test and six neighboring tests passed
in the assistant's clean GNU Fortran 14.2.0 Debug build. CMake was 3.31.6.
No user's workstation or running simulation was accessed.

```
cmake -S . -B build-review -G Ninja -DCMAKE_BUILD_TYPE=Debug -DPELEF_ENABLE_TESTS=ON
cmake --build build-review --parallel 4 --target test_nasa7_thermo test_mixture_thermo test_amr_reactive_3d_checkpoint test_reactive_eb_3d_checkpoint test_full_h2o2_jacobian test_implicit_h2o2_reactor test_pressure_dependent_kinetics
ctest --test-dir build-review --no-tests=error --output-on-failure -j 2 --tests-regex '^unit_(nasa7_thermo|mixture_thermo|amr_reactive_3d_checkpoint|reactive_eb_3d_checkpoint|full_h2o2_jacobian|implicit_h2o2_reactor|pressure_dependent_kinetics)$'
```

Result: **7/7 passed, zero failures and zero skips**. The three new rejection
scenarios are inside the existing NASA7 test and are not counted as three
additional CTest registrations. Project-contract and whitespace checks passed.
The exact tested blobs subsequently published are:

- `src/physics/nasa7_thermo_mod.F90`: `7526898120b3b10d75ff5fe208599cd81fe33395`
- `tests/unit/test_nasa7_thermo.F90`: `1c177bcab06c54d3f203dd22d02a4379df2f9c19`

## Qualification boundary

These seven focused tests are not a full matrix or physical validation.
The existing complete serial, MPI, selected CVODE/Cantera and tests-disabled
installed-program workflows must qualify the repaired source. The authoritative
completed run identities and results are recorded in PR #102. A successful
run on the parent source does not qualify this failure-status correction.
