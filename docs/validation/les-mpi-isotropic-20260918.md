# MPI static LES isotropic stress — 2026-09-18

Implementation `b2ba2e58c08c49a35c80780633e0e5960510d388`, based on `81e2ea0`.

`advance_mpi_spray_les` accepts trailing optional `ci`, default zero. It validates
finite CI in `[0,1]` collectively and checks exact rank agreement before halo
communication or candidate updates. Both periodic and physical-boundary paths
forward it to the shared ghosted RHS at both SSPRK2 stages. The existing
transactional state/temperature commit and boundary ledger remain in effect.
Coupled-spray/application controls are still pending. Reference formulas and
interval-estimate limitations are in [the serial record](les-isotropic-20260918.md).

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug, existing configurations unchanged.
Fixed tests-enabled `test_les_mpi_models` and selected H2O2 tests-disabled
`pelef_selected_spray_runtime` built successfully. Selected runtime was not run.
Focused CTest **2/2 passed**, 2.90 seconds, zero failures/skips:
`spray_les_mpi_models_2` (1.50 s), `spray_les_mpi_models_4` (1.35 s).

Each run exercises Smagorinsky/WALE/Vreman with CI zero, CI 0.1 plus eddy
coefficient 0.17, and CI 0.1 with eddy coefficient zero, for periodic and
outflow/isothermal moving-wall boundaries. State relative differences stay
below `2e-12` with unit floor, temperature absolute differences below `1e-8 K`,
and boundary changes within `1e-10` relative with `1e-12` scale floor.
These are existing parity tolerances, unchanged. Checks also cover nonzero
isotropic evolution, rank-dependent CI rejection at zero dt and zero eddy
coefficient, one-rank NaN CI, excessive dt, and model disagreement/invalidity;
all rejected updates preserve state and temperature exactly.

Provenance, Python syntax, project contract and whitespace checks passed.
No builds/tests failed. No full matrix, external solver parity or physical
validation was performed. Local slabs still require at least two x cells.

## Commands and artifacts

Logs: `build/les-mpi-isotropic-checks/{provenance,build,build-selected,ctest}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_models --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_mpi_models_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
