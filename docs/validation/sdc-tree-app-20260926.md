# Static serial AMR SDC application — 2026-09-26

Implementation: `8f11c535` (base `140481dc`). PeleF 0.246.0.
Pinned PeleC: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
GNU Fortran 13.3.0, Ninja Debug; existing fixed tests-enabled and selected
tests-disabled build directories. No compiler/backend configuration change.

Fixed and selected frontends now select the common-step tree SDC integrator.
The driver keeps reaction estimates across accepted steps and separately sums
external-force inventories. Failed attempts retain the integrator's transactional
state. Split output remains unchanged. SDC restart, regridding and recursive
subcycling are explicitly rejected until their state contracts are implemented.

## Completed focused checks

Artifacts: `build/sdc-tree-app-checks/` (local, not committed).

- Mechanism provenance preflight passed: `h2o2_full.json`, bundle SHA
  `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Built `pelef_amr_tree_reactive_3d` in `build/spray-amr-reactive-debug` and
  `pelef_amr_tree_reactive_3d_selected` in `build/spray-amr-reactive-selected-debug`:
  `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build BUILD --target TARGET --parallel 2`.
  Logs: `build.log`, `build-selected.log`; final fixed rebuild: `build-final.log`.
- Two application runs passed, each with two levels, two steps to 1e-9 s,
  two SDC iterations, reactions, all configured transport and nonzero forcing.
  Inputs/logs: `input.nml`, `selected.nml`, `topology.dat`, `run.log`, `run-selected.log`.
  Command: `BUILD/TARGET INPUT build/sdc-tree-app-checks/topology.dat` with
  `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset. Fixed fixture is committed as
  `cases/amr_tree_reactive_3d/sdc.nml` and `sdc.dat` (output path changed).
- Fixed/selected state and budget CSVs are byte-identical. Maximum relative
  Euler inventory residual after boundary and external-source subtraction:
  **1.5088880919823904e-16**. See `comparison.txt`. Transport theta: 1.
- Three unsupported-input checks passed: subcycling, checkpoint writing,
  automatic regridding each failed with its specific SDC diagnostic.
  Corresponding `.nml` and `.log` files are retained.
- Python syntax for `tools/check_project_contract.py`, project contract and
  `git diff --check` passed.

Build warnings remain. No CTest suite, physical-boundary AMR case, chemistry-off
application case, physical validation, installation audit or full solver parity
was run in this focused change. The successful numerical runs preceded a final
whitespace-only wrapping of the split call; the fixed executable was rebuilt.
