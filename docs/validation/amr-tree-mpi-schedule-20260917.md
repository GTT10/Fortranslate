# Scheduled MPI tree application regridding — 2026-09-17

## Source and scope

Base: `c2b6b64e5980c790ec217f7dd40a836c19446491`, local `main`.
Fixed/selected MPI frontends now accept the serial schedule as a third argument.
The driver connects sparse remapping, deterministic compensated inventories,
exact event clipping, geometry/owner replacement and schedule-aware restart.
No mechanism source, numerical kernel, tolerance or golden hash changed.

Final source SHA-256 (driver comment and checker wording cleaned after compilation):

- `src/driver/mpi_amr_patch_tree_application_3d_mod.F90`:
  `9d65b114826fa7d47f5a58729bc4182d6c15e104849af9faabc7b38c8bcaa9cf`
- `app/pelef_mpi_amr_tree_reactive_3d.F90`:
  `420b3334f4d4ba446d49a7d01a2ec3946ce305151e07e77f4c8f6b36a9806cda`
- `app/pelef_mpi_amr_tree_reactive_3d_selected.F90.in`:
  `9a96167ad2b2a1a1b3f2c5362ed462e7c1a8b1d0665e5359b960c1755852aed2`
- `tools/check_mpi_amr_tree_restart_3d.py`:
  `899676b80915e56b7745a00288018f2b089a810a9db81a517821beebd5c84b1c`

## Verification

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, existing clean system-MPI configuration
`build/selected-tree-mpi-debug`; compiler/backend settings unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 4 --target pelef_mpi_amr_tree_reactive_3d test_selected_elementary_mpi_amr_tree_3d test_selected_elementary_amr_tree_3d pelef_amr_tree_reactive_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_mpi_amr_tree_(schedule_3d|restart_3d)$' --output-on-failure
python3 -m py_compile tools/check_mpi_amr_tree_restart_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The first build rejected an overlong selected-frontend call; splitting the source
line corrected it. Rebuild and MPI ABI 40 audits passed. Existing warnings remain.
Final CTest **2/2 passed**, 107.17 s total: static restart 55.25 s, scheduled
regrid/restart 51.89 s; no failures/skips. Python syntax, provenance, project
contract and whitespace checks passed. Fixed/selected scheduled runs report
maximum normalized remap conservation error 1.5696e-16 (limit 5e-12).

The scheduled fixture changes to ratio-3 branches at 2.5e-10 s, removes all children
at 5e-10 s, and restores the original branched tree at 1e-9 s. It performs three
steps/events. For fixed and selected elementary mechanisms the checker compares:

- Serial versus 2-rank uninterrupted CSV/budget bytes.
- 2-rank checkpoint after event 1 versus 3-rank continuation, with identical output.
- MPI checkpoint/serial continuation and serial checkpoint/1-rank continuation;
  checkpoint bytes match as well as final outputs.
- Stop after event 2 (root-only), then resume with 3 ranks, including two empty
  owners, through the final refinement event.
- Changed schedule/physics, truncated checkpoint, inconsistent finite ledger,
  output alias and malformed namelist rejection, preserving final outputs.

The existing static restart regression checks the same continuation controls
without schedule events. These are focused application comparisons, not a full
configuration matrix. No installed-prefix audit, scaling study or external/physical
validation was run.

Logs: `build/amr-tree-mpi-schedule-checks/{build,ctest}.log`.
Artifacts: `build/selected-tree-mpi-debug/tests/mpi_amr_tree_schedule_3d/{fixed,selected}/`
and the static `mpi_amr_tree_restart_3d/` sibling directory.

## Limits

Automatic temperature-driven MPI adaptation and subcycling remain open. Schedules
are replicated input metadata. Remap payloads are sparse; the inventory check
broadcasts sum/correction vectors per node. Root gathers remain at initialization,
checkpoint and final output. The existing path and checkpoint-I/O limitations apply.
