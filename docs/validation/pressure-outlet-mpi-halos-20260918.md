# Sparse MPI pressure outlet halos — 2026-09-18

Implementation `bcfc5e5`, base `57a90ce`.

`build_mpi_amr_sparse_boundary_halos_3d` accepts optional species thermodynamics
and `outlet_pressure(6)`. Pressure input requires species; finite targets in
[0,1e15] may be positive only on outflow faces. Optional-presence flags and
pressure values are checked collectively before EOS processing.

Transverse pressure ghosts are applied around exchanged neighbor columns;
x-edge pressure ghosts use the nearest two real global columns. This preserves
one-column-per-rank support without a full-level gather. Face-order corner
composition matches the serial helper. Candidate halos publish only after
collective success; rejected calls preserve caller outputs.

This connects the halo API only. MPI hydro/transport callers, fine-patch
pressure interpolation, application input and patch trees remain pending.

## Focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6, existing Debug/Ninja fixed tests-enabled
and selected H2O2 tests-disabled builds. Added cases compare every local halo
entry exactly with the corresponding serial full-grid slice, including mixed
walls/inflow/outflow and pressure corners. A rank-local invalid pressure
assignment must reject collectively without changing output arrays. Existing
physical metadata rejection checks are included. The `--halos-only` mode exits
after initialization and physical halo checks; full sparse AMR tests remain
available in the default mode.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_reactive_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/mpiexec -n 2 build/spray-amr-reactive-debug/test_mpi_amr_reactive_3d --halos-only
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/mpiexec -n 4 build/spray-amr-reactive-debug/test_mpi_amr_reactive_3d --halos-only
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_reactive_3d_selected --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/pressure-outlet-mpi-halo-checks/`. This is focused halo
parity/rollback evidence, not an MPI pressure-outlet flow or physical validation.

Final focused launches **2/2 passed** (2 and 4 ranks, exit 0), zero skipped
launches. Fixed and selected executable builds, selected MPI ABI audit,
provenance and static checks passed. The original full 2-rank launch and the
first focused-mode attempt (whose branch was initially placed after expensive
unrelated tests) were stopped to limit verification; neither counts as passed.
The final focused branch runs immediately after hierarchy initialization.
No tolerance changes. Runtime logs retain X authorization messages and
underflow/denormal notes; both final processes exited successfully.

Final logs: `build-focused.log`, `build-selected-final.log`, `np2-focused.log`,
`np4-focused.log`. Earlier build/run logs are retained in the same directory.
