# Sparse MPI pressure outlet Strang split — 2026-09-18

Implementation `3980e75`, base `66cb298`.

`advance_mpi_amr_sparse_strang_3d` accepts optional `outlet_pressure(6)`,
validates target values and includes them in the collective contract before
chemistry. The private candidate then advances reaction/hydro/reaction with
pressure passed into hierarchy hydro. Rejected stages retain the existing
state and diagnostic rollback. MPI molecular transport, full transport split
and application pressure input remain pending.

## Focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6, existing Debug/Ninja fixed tests-enabled
and selected H2O2 tests-disabled builds.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_reactive_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/mpiexec -n 2 build/spray-amr-reactive-debug/test_mpi_amr_reactive_3d --pressure-strang-only
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/mpiexec -n 4 build/spray-amr-reactive-debug/test_mpi_amr_reactive_3d --pressure-strang-only
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_reactive_3d_selected --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The focused fixture uses an upper-corner patch with 102325 Pa upper-face
targets, active chemistry, PCM/characteristic PLM and a 1e-10 s interval.
Gathered coarse/fine states and temperatures and reflux diagnostics must equal
the serial Strang reference exactly. Existing surface, Euler and elemental
checks retain their tolerances. A rank-local invalid inlet rejects the private
candidate while preserving both levels and zeroing boundary/reflux diagnostics.
Selected runtime is compiled only. This does not establish ignition or acoustic
validation.

Artifacts: `build/pressure-outlet-mpi-strang-checks/` contains `provenance.log`,
`build.log`, `build-selected.log`, `np2.log`, and `np4.log`.

Final focused launches **2/2 passed**, 2 and 4 ranks, exit 0; no skipped,
failed or interrupted launches. Fixed/selected builds, selected MPI ABI audit,
provenance and static checks passed. No tolerance changes. Logs retain X
authorization messages and underflow/denormal notes despite successful exits.
