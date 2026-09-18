# MPI two-level pressure outlet application — 2026-09-18

Regrid implementation `8e99bd9`, application/tests `46f3b5b`, base `74e27aa`.

Fixed and selected MPI two-level AMR applications enable pressure input,
forward targets into complete updates and provide them during distributed
regridding. Regrid validates target agreement before remapping and uses
pressure-aware coarse halos for interpolation. The existing shared checkpoint
boundary identity already stores all six targets, so the MPI path also rejects
changed or removed targets. Patch-tree pressure input remains pending.

## Focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6, existing Debug/Ninja fixed tests-enabled
and selected H2O2 tests-disabled configurations. Selected MPI executable built;
selected runtime not executed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_mpi_amr_reactive_3d pelef_amr_reactive_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_reactive_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^regression_mpi_amr_pressure_outlet_restart_3d$' --output-on-failure
python3 -m py_compile tools/check_amr_physical_boundary_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The existing pressure fixture uses inert hydro and molecular transport, moving
isothermal x-low wall, 102325 Pa x-high target and temperature-triggered
boundary-touching regridding. A serial full reference is compared against a
4-rank checkpoint followed by a 2-rank restart: coarse/fine CSVs must match
byte-for-byte, and surface ledger differences must meet the existing MPI
relative 2e-11 plus 1e-25 gate with common momentum scale. Changed pressure,
removed pressure and changed wall temperature must reject without output CSVs.
Six executable launches. This is application/restart equivalence, not physical
validation or scaling evidence.

Artifacts: `build/pressure-outlet-mpi-app-checks/` contains provenance,
fixed/selected build logs and `ctest.log`.

Final CTest **1/1 passed**, zero skips, 69.65 s total. Fixed/selected builds,
selected MPI ABI audit, provenance, Python syntax, project contract and diff
checks passed. No failed/interrupted runs or tolerance changes.
Final fixture directory:
`build/spray-amr-reactive-debug/tests/mpi_amr_pressure_outlet_restart_3d/run-8cwu06li/`.
