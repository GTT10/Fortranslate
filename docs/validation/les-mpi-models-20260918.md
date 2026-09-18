# MPI LES model dispatch — 2026-09-18

Implementation `a5f1855`, base `19e62f6`.

`advance_mpi_spray_les` accepts optional trailing `model` with `smagorinsky`,
`wale`, or `vreman`. Omission resolves to Smagorinsky. The normalized model ID
participates in collective configuration agreement before halo communication
or candidate evolution. Invalid names reject collectively, including zero-dt/
zero-coefficient calls. Both periodic and physical-halo branches forward the
selected model to each RK2 stage. Existing diffusion bounds, global boundary
ledger and state/temperature publication are retained. Coupled spray and CLI/
checkpoint model selection still require connection.

## Focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6, `/usr/bin/mpiexec`; reused Debug/Ninja fixed
tests-enabled and selected H2O2 tests-disabled configurations.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_models --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_mpi_models_(2|4)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The dedicated 8x4x4 crossed-velocity test covers all three models at two/four
ranks, including two-cell local x extents. Periodic and outflow/moving-isothermal
y-wall cases compare conserved states (relative 2e-12, unit floor), temperature
(absolute 1e-8) and boundary inventories (relative 1e-10, floor 1e-12) with serial.
A last-rank invalid name and a last-rank different valid model reject while
preserving every local state/temperature entry. These gates follow the existing
MPI physical SGS comparisons; no existing tolerances changed.

CTest **2/2 passed**, zero skips, 1.58 s total (0.78/0.77 s for two/four ranks).
Fixed/selected builds, provenance, Python syntax, project contract and diff checks
passed. No failed or interrupted runs. Selected compilation is separate from
fixed-mechanism runtime parity; no physical-validation or scaling claim.

Logs: `build/les-mpi-model-checks/{provenance,build,build-selected,ctest}.log`.
