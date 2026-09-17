# MPI single-depth hydro SSPRK2 — 2026-09-17

## Source and scope

Base: `aabeec1f3e0defd63664a2e5c4a125395d4b4052`, local `main`.
Adds a separate MPI hydro level-step API, fixed/selected runtime wiring and focused
serial parity tests. Existing synchronous hydro and serial numerical kernels are
unchanged. The caller retains responsibility for level scheduling and reflux.

Verified SHA-256:

- `src/parallel/mpi_amr_patch_tree_hydro_level_3d_mod.F90`:
  `ec5ba51c85d9a9db9a964442f5ddad0fd097a0a6fd7dc45b87f991c5fb0f3c63`
- `tests/unit/test_mpi_amr_patch_tree_hydro_level_3d.F90`:
  `90cff75f624acfc0c5dcb1c2962f677ba678b7e3ed60ddcf736d6b397056d9cf`

## Completed checks

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, existing clean system-MPI configuration
`build/selected-tree-mpi-debug`; compiler/backend options unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 4 --target pelef_mpi_support pelef_selected_mpi_amr_tree_reactive_3d_runtime
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 4 --target test_mpi_amr_patch_tree_hydro_level_3d pelef_selected_mpi_amr_tree_reactive_3d_runtime
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_mpi_amr_patch_tree_hydro_level_3d_[23]$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Both builds passed. CTest **2/2 passed**, 9.35 s total: 2 ranks 4.69 s, 3 ranks
4.64 s; no failures/skips. Provenance, project contract and whitespace checks
passed. No Python sources changed. Existing compiler warnings remain.

The three-level H2/O2 fixture has ratios 2 then 3 and tests depths 0/1/2 with
PCM/characteristic PLM, periodic, outflow, slip walls and reservoir inflow. Owner
maps include empty ranks and all payload on one rank. Coarse endpoint state
changes by 2%; stage fractions are [0.1,0.25] then [0.6,0.75], truncated to the
active depth. Endpoint fields at/above active depth are absent.

Every owned conserved/temperature field matches serial exactly; inactive fields
remain equal to their inputs. All three integrated face-flux arrays match serial
exactly on active owners, while inactive flux slots are absent. Invalid solver
and rank-dependent second-stage fractions are collectively rejected with prior
state, temperature and integrated flux preserved. The latter checks rollback
after the first stage has completed privately. No tolerances or hashes changed.

Build/test logs: `build/amr-tree-mpi-hydro-level-checks/{build,ctest}.log` (the build
log contains the final incremental test build; the runtime build completed first).
Binary: `build/selected-tree-mpi-debug/tests/test_mpi_amr_patch_tree_hydro_level_3d`.

## Limits

This API neither averages down nor synchronizes/refluxes face fluxes. It does not
compute/enforce CFL or retry dt. Species coefficient identity is a caller contract.
MPI transport level steps and recursive integration/application subcycling remain
open. No full matrix, installed-prefix audit, scaling study, external solver parity
or physical validation was run.
