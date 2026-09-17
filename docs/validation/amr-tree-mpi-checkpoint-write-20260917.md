# Owner-local MPI tree checkpoint saving — 2026-09-17

## Source and implementation

Base `f7d66388315749488da63359392070274f12fc5d`, local `main`.
Adds MPI checkpoint saving with owner-local EOS/temperature and parent-average
preflight, collective metadata/context agreement and checked byte offsets.
Each owner writes its state/temperature records directly with MPI-IO using
bounded chunks. Rank zero writes the header/footer. Fixed/selected application
save dispatch no longer gathers fields; serial V1 format and restart identity
are preserved. Loading remains the root serial reader plus scatter.
No serial numerical algorithm, mechanism source, tolerance or golden hash changed.

SHA-256:

- `src/parallel/mpi_amr_patch_tree_checkpoint_3d_mod.F90`:
  `7f49fb7c6729e33aa0ed9cfd0a7ba9527fbaf1db6a37ab6e90dfe49822b70c68`
- `src/driver/mpi_amr_patch_tree_application_3d_mod.F90`:
  `448cafe0409e3d4c7436ce0c46b88207a6cbd629594398907d962aad9e84d73d`
- `tests/unit/test_mpi_amr_patch_tree_checkpoint_3d.F90`:
  `362cadd0440c59c7323d3249cba3b8f80fcbe84c1ae40f183fe25e554fd705dc`

## Commands

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, existing
`build/selected-tree-mpi-debug`; no compiler/backend changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_mpi_amr_patch_tree_checkpoint_3d pelef_mpi_amr_tree_reactive_3d test_selected_elementary_mpi_amr_tree_3d -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(unit_mpi_amr_patch_tree_checkpoint_3d_[23]|regression_mpi_amr_tree_(restart|weighted_adaptation)_3d)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

## Completed results

Fixed/selected application and unit builds passed. CTest **4/4 passed**,
98.62 s total: static restart 55.41 s, weighted adaptation 42.26 s, unit
2/3 ranks 0.45/0.46 s. No failures/skips. Provenance, project contract and
whitespace checks passed. Existing compiler warnings remain.

## Coverage

The unit fixture has parent/child owners on different ranks and an empty owner
at three ranks. Incorrect child temperature, a covered-parent state inconsistent
with its children, invalid time/step metadata and rank-dependent context must
reject collectively while preserving prior file content and length. A valid
owner-written checkpoint must load through the serial reader with exact U/T.

Existing fixed/selected application tests compare checkpoint bytes with serial,
continuous/restarted CSV and budgets, 2-to-3-rank and serial interchange. They
cover common-step reacting/transport execution and recursive PLM automatic
adaptation with weighted ownership, including a policy change on restart.
Existing corrupted-file/context/ledger guards remain enforced.

Logs: `build/amr-tree-mpi-checkpoint-write-checks/{build,ctest}.log`.
Artifacts: `build/selected-tree-mpi-debug/tests/mpi_amr_tree_{restart,weighted_adaptation}_3d/`
and `mpi_tree_checkpoint_preflight_{2,3}.chk` in the test directory.
No Python sources changed.

## Limits

Distributed validation makes temporary owner-local field copies. Header metadata
is replicated. Actual I/O failure after open/truncation can leave a partial file;
only preflight failure preserves the destination. Checkpoint loading still holds
all fields on rank zero. No filesystem throughput/scaling benchmark, full matrix,
installed audit, external-solver parity or physical validation was run.
