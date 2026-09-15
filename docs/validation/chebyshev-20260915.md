# Chebyshev rates — 2026-09-15

Source: parent `90f2a5a` plus this commit's runtime, importer/generator,
checkpoint/MPI metadata, tests and documentation. GNU Fortran 13.3 Debug,
`build/port-spray-debug`; Cantera 3.2.0 oracle.

## Implemented scope

Reaction kind 10 implements the temperature/pressure Chebyshev expansion in
[Cantera 3.2 rate parameterizations](https://cantera.org/3.2/reference/kinetics/rate-constants.html#chebyshev-reaction-rate-expressions).
The coefficient matrix has at most eight temperature and eight pressure terms.
The reduced coordinates use inverse temperature and log10 pressure; physical
pressure is R*T*sum(C). Polynomial recurrence supplies both the rate and its
analytic concentration derivative. All species have the same effective-rate
pressure derivative; the existing mass-action terms remain separate.

Positive ordered fit bounds and finite rectangular coefficients are required.
Zero total concentration is rejected. Evaluation extrapolates outside the fit
range, with the existing absolute natural-log rate limit of 700; this does not
validate the fit outside its specified range. Extremely ill-conditioned fits
and extreme extrapolation are not qualified by the numerical fixture.

The importer reads normalized coefficients from Cantera's rate API, avoiding
local input-unit metadata ambiguity. The generator emits explicit coefficients
and bounds. AMR/EB contexts serialize shape, four bounds and the fixed matrix
only for the new reaction kind. Dense/sparse MPI comparisons include all fields.
The additional fixed storage is 68 reals and two integers per reaction; MPI
comparison adds 70 reals. Large-mechanism memory/scaling is not benchmarked.

## Focused evidence

- Python portability: 13 tests passed, including Chebyshev API coefficient
  selection, generation, invalid bounds, ragged/empty/oversized tables and NaN.
- Generator regression: 12 tests passed in 22.567 s; old generated output remains.
- Fortran kinetics and AMR checkpoint variants: 8 tests passed in 0.65 s.
  The new variant round-trips and rejects a changed mixed temperature/pressure
  coefficient transactionally.
- The oracle constructs synthetic H+OH<=>H2O with a 3 by 3 matrix, fit ranges
  300--2500 K and 1e3--1e7 Pa, exports cm/mol/atm YAML, imports and generates a
  Fortran loader, then compiles a bounds/FPE-checked rate probe. Forty-two states
  span 250, 300, 700, 1600, 2500, 3000 K and 1e2, 1e3, 3e4, 1e5, 7e6, 1e7,
  1e8 Pa. Interior, endpoint and extrapolated rates have maximum relative error
  1.1768364061e-14. Every concentration derivative agrees with fourth-order
  Cantera finite differences within 3.7133629505e-11 relative error.
- The compiled probe also checks a constant 1 by 1 polynomial and rejection of
  zero pressure, an oversized matrix and a degenerate temperature interval.
- GNU/OpenMPI `pelef_mpi_amr_reactive_3d` rebuilt successfully in the existing
  clean `build/port-regrid-mpi`; new Chebyshev MPI execution and EB restart were
  not run. This is rate/checkpoint evidence, not reactor or physical validation.
- Provenance preflight, project contract, Python syntax and whitespace passed.

Artifacts: `build/chebyshev-check/chebyshev-f6bi94r2/` contains YAML/JSON,
generated Fortran, probe, rate output and `result.json`. Source SHA256:
`c957a1c055f48756881e7431d94bb99557cb4927d4c86db13d7e368fafa4d9c2`.
Build logs: `build/port-spray-debug/chebyshev-build.log` and
`build/port-regrid-mpi/chebyshev-build.log` (existing warnings, no build errors).

```sh
cmake --build build/port-spray-debug --target test_elementary_kinetics test_amr_reactive_3d_checkpoint --parallel 8
ctest --test-dir build/port-spray-debug --output-on-failure -R '^unit_(elementary_kinetics|amr_reactive_3d.*checkpoint)$'
python3 -m unittest discover -s tests/python -p test_generate_mechanism_bundle.py
python3 -m unittest discover -s tests/python -p test_mechanism_portability.py
build/cantera-venv/bin/python tools/check_chebyshev_rates.py --build-dir build/port-spray-debug --work-dir build/chebyshev-check
cmake --build build/port-regrid-mpi --target pelef_mpi_amr_reactive_3d --parallel 8
python3 tools/check_mechanism_provenance.py --bundle build/chebyshev-check/chebyshev-f6bi94r2/chebyshev.json --source build/chebyshev-check/chebyshev-f6bi94r2/chebyshev.yaml
python3 tools/check_project_contract.py --source .
python3 -m py_compile tools/generate_elementary_mechanism.py tools/ingest_cantera_mechanism.py tools/check_chebyshev_rates.py tests/python/test_mechanism_portability.py
git diff --check
```
