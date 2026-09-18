# Sparse MPI pressure outlet coarse hydro — 2026-09-18

Implementation `94e6c9a`, base `be03ab5`.

`advance_mpi_amr_sparse_boundary_hydro_3d` forwards optional
`outlet_pressure(6)` and species thermodynamics into both RK2 halo builds.
The existing collective metadata validation and private-candidate publication
handle rank-local rejection. Surface accounting continues to use accepted
stage fluxes. This connects coarse hydro only; fine-patch pressure ghosts,
full hierarchy dispatch, transport and MPI application input remain pending.

## Completed focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6, existing Debug/Ninja fixed tests-enabled
and selected H2O2 tests-disabled configurations.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_reactive_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/mpiexec -n 2 build/spray-amr-reactive-debug/test_mpi_amr_reactive_3d --coarse-hydro-only
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/mpiexec -n 4 build/spray-amr-reactive-debug/test_mpi_amr_reactive_3d --coarse-hydro-only
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_mpi_amr_reactive_3d_selected --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Focused launches **2/2 passed**, exit 0, no skipped or interrupted launches.
Both builds, selected MPI ABI audit, provenance and static checks passed.
The fixture covers PCM and characteristic PLM, zero and 102325 Pa targets,
and mixed inflow/wall/outflow faces. Every local state, temperature and tested
mean face flux equals its serial reference exactly. Global surface inventory
matches at the existing relative 2e-11 tolerance with floor 1e-20. A rank-local
invalid inlet with pressure enabled preserves state/temperature and leaves a
zero boundary ledger. No tolerance changes.

The focused mode runs after hierarchy initialization and skips unrelated
full-suite work. Selected runtime was compiled, not executed. Logs retain
X authorization messages and underflow/denormal notes; successful exit status
is recorded independently. No acoustic-reflection or physical-validation claim.

Artifacts: `build/pressure-outlet-mpi-hydro-checks/` contains `provenance.log`,
`build.log`, `build-selected.log`, `np2.log`, and `np4.log`.
