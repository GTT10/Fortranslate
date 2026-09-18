# Sparse MPI pressure outlet transport/full split — 2026-09-18

Implementation `5ea2b2c`, base `3f38e77`.

Sparse MPI transport and full reaction/transport/hydro wrappers now accept
`outlet_pressure(6)`, validate it and include targets in their collective
contracts. Targets reach coarse physical transport fluxes, coarse start/end
interpolation halos, fine transport ghosts, hydro and the transport-disabled
Strang dispatch. Existing limiting, subcycling, flux accounting and rollback
remain intact. MPI application input/restart/regridding and patch trees are
still pending.

## Focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6, existing Debug/Ninja fixed tests-enabled
and selected H2O2 tests-disabled configurations.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_reactive_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/mpiexec -n 2 build/spray-amr-reactive-debug/test_mpi_amr_reactive_3d --pressure-transport-only
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/mpiexec -n 4 build/spray-amr-reactive-debug/test_mpi_amr_reactive_3d --pressure-transport-only
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_reactive_3d_selected --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The focused upper-corner fixture uses transport-only and reacting full updates
at 1e-10 s, y-high/z-high targets of 102325 Pa and moving isothermal x walls.
Gathered coarse/fine states and temperatures and reflux diagnostics must match
serial results exactly; existing surface, component/Euler and elemental checks
retain their tolerances. Invalid inlet metadata rejects collectively and
preserves both hierarchy levels, boundary ledger, reflux and limiter diagnostics.
Selected runtime is built only. This is not physical or acoustic validation.

Artifacts: `build/pressure-outlet-mpi-transport-checks/` contains `provenance.log`,
`build.log`, `build-selected.log`, `np2.log`, and `np4.log`.

Final focused launches **2/2 passed**, 2 and 4 ranks, exit 0; no skipped,
failed or interrupted launches. Fixed/selected builds, selected MPI ABI audit,
provenance, Python syntax, project contract and diff checks passed. No tolerance
changes. Runtime logs retain X authorization and underflow/denormal messages
while exiting successfully.
