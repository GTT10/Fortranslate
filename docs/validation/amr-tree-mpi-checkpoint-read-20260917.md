# Indexed owner-local MPI checkpoint loading — 2026-09-17

## Source and implementation

Base `bfe3f5178dabe522c5350113ce921b5531d0faac`, local `main`.
The serial header parser is shared with the MPI reader. Rank zero streams through
V1 cell records using one-cell scratch storage, builds node-position offsets and
checks footer/EOF. Owners then seek to their ranges and read only their assigned
U/T payloads. Collective EOS/parent-average preflight is shared with saving.
The application checks geometry/counters and distributed inventories before
publishing candidates; the full root snapshot, root field scatter and separate
adaptive-geometry broadcast are removed. Rank/policy-neutral V1 compatibility,
original numerical kernels, mechanisms and validation tolerances are retained.

SHA-256:

- `src/io/amr_patch_tree_checkpoint_3d_mod.F90`:
  `168f9cdae45e5bd36f11b85d104526179edabfd0ed507477e9f72a0f53d066f5`
- `src/parallel/mpi_amr_patch_tree_checkpoint_3d_mod.F90`:
  `0d85b06ed1220d99550d48670547a5531da6d76bf1aac3e27fdf42f7c58d12e0`
- `src/driver/mpi_amr_patch_tree_application_3d_mod.F90`:
  `6ad39c850286aa6a6e56dff105c841dd29a55efd9dcb60fb22584a77f07dff82`
- `tests/unit/test_mpi_amr_patch_tree_checkpoint_3d.F90`:
  `98fa4a59729847083a160e064db9818f114ee8f8b043c793b3d93a174d5b2593`

## Completed checks

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, existing
`build/selected-tree-mpi-debug`; compiler/backend settings unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target pelef_mpi_amr_tree_reactive_3d test_selected_elementary_mpi_amr_tree_3d pelef_amr_tree_reactive_3d test_selected_elementary_amr_tree_3d test_mpi_amr_patch_tree_checkpoint_3d -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(unit_mpi_amr_patch_tree_checkpoint_3d_[23]|regression_mpi_amr_tree_(restart|weighted_adaptation)_3d)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed/selected serial and MPI builds passed. CTest **4/4 passed**, 100.03 s total:
static restart 56.12 s, weighted adaptation 42.96 s, checkpoint 2/3 ranks
0.45/0.46 s. No failures/skips. Provenance, project contract and whitespace
checks passed. No Python sources changed; existing compiler warnings remain.

Application checks compare exact serial/MPI checkpoint bytes and continuous/
restarted CSV/budget bytes, including 2-to-3-rank and serial interchange, full
common-step reacting/transport and recursive PLM hydro-only automatic adaptation.
Weighted-to-round-robin restart remains exact. Existing context/geometry/counter,
truncation and boundary-ledger guards remain enforced.

The checkpoint unit adds changed ownership on read, owner-only allocation and
exact U/T assertions. It corrupts only a temperature record through unformatted
positioned output, explicitly checking that file size is unchanged; distributed
EOS validation must reject the candidate while preserving prior owner vector and
U/T. Existing save preflight and serial-reader checks also run. The refined
corruption fixture was rebuilt before these unit tests executed.

Logs: `build/amr-tree-mpi-checkpoint-read-checks/{build,test-build,ctest}.log`.
Artifacts: `build/selected-tree-mpi-debug/tests/mpi_amr_tree_{restart,weighted_adaptation}_3d/`
and `mpi_tree_checkpoint_preflight_{2,3}.chk` in the test directory.

## Limits

Rank zero still performs a serial indexing scan and payloads are read again by
owners. This removes global-state memory and field scattering, not the serial
scan or duplicated read traffic. Metadata/index arrays are replicated and
validation uses temporary owner-local copies. No indexed format, parallel I/O
throughput benchmark, full matrix, installed audit, external-solver parity or
physical validation was added. Files must remain unchanged during a read.
