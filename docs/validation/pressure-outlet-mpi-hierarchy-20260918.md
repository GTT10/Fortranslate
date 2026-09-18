# Sparse MPI two-level pressure outlet hydro — 2026-09-18

Implementation `d17e216`, base `d667ae3`.

`advance_mpi_amr_sparse_hydro_3d` accepts `outlet_pressure(6)` and includes
normalized target values in its collective contract. It forwards targets to
coarse hydro, coarse start/end interpolation halos and both fine RK2 stages.
The fine extended-state builder applies pressure only at physical root faces
touched by the patch, using local or exchanged adjacent fine donors and global
fine coordinates. Face order matches the serial corner convention.

Subcycling, accepted-flux accumulation, reflux, average-down and composite
surface accounting retain their existing algorithms and transactional commit.
Strang/full reaction wrappers, molecular transport and MPI application input
are still pending; this increment covers hierarchy hydro.

## Focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6, existing Debug/Ninja fixed tests-enabled
and selected H2O2 tests-disabled configurations.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_reactive_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/mpiexec -n 2 build/spray-amr-reactive-debug/test_mpi_amr_reactive_3d --pressure-hierarchy-only
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/mpiexec -n 4 build/spray-amr-reactive-debug/test_mpi_amr_reactive_3d --pressure-hierarchy-only
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_reactive_3d_selected --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The focused mode runs PCM/characteristic PLM on interior, low-corner,
high-corner and fully refined patches, with targets on all three upper faces.
Gathered coarse/fine states and temperatures and the maximum reflux diagnostic
must exactly match serial two-level hydro. Surface ledgers use the existing
relative 2e-11 tolerance with floor 1e-20. Composite conservation and rank-local
invalid-inlet rollback checks remain. Selected runtime is compiled only.

Artifacts: `build/pressure-outlet-mpi-hierarchy-checks/`, containing provenance,
initial/final fixed and selected build logs, `np2.log`, and `np4.log`.
This is focused numerical equivalence, not acoustic or physical validation.

Final focused launches **2/2 passed**, exit 0 at 2 and 4 ranks; no skipped,
failed or interrupted launches. Fixed/selected builds, selected MPI ABI audit,
provenance, Python syntax, project contract and diff checks passed. Runtime
logs retain X authorization messages and underflow/denormal notes; these did
not prevent successful exits. No tolerance changes.
