# Static MPI tree checkpoint/restart — 2026-09-17

## Source and scope

Base: `c7371be5b6f6695d30658e72929ceab4b59b145a`, local `main`.
Connects fixed/selected MPI tree checkpoint controls to the existing serial
checkpoint schema and validator. Root-only checkpoint I/O, owner redistribution
and metadata broadcasts support changed-rank continuation. No numerical kernels,
mechanism data, tolerances or golden hashes changed.

Verified SHA-256:

- `src/driver/mpi_amr_patch_tree_application_3d_mod.F90`:
  `734434258814abe94648ca270f59ec2aef3364d5111f6d69abeba699a86b9601`
- `app/pelef_mpi_amr_tree_reactive_3d_selected.F90.in`:
  `4101a2f75cca13df2d187bc96882b1fb6e7b9c6d393995de056628ac6b9a47f4`
- `tools/check_mpi_amr_tree_restart_3d.py`:
  `a0a0be5ec3068e0e147df3d02ae40f1700a1990b00e00929a576f4bd3766dd56`

## Completed verification

GNU Fortran 13.3.0 Debug, OpenMPI 4.1.6, existing clean system-MPI configuration
`build/selected-tree-mpi-debug`; compiler/backend settings unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 4 --target pelef_mpi_amr_tree_reactive_3d test_selected_elementary_mpi_amr_tree_3d test_selected_elementary_amr_tree_3d pelef_amr_tree_reactive_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_mpi_amr_tree_(restart_3d|selected_application_3d)$' --output-on-failure
python3 -m py_compile tools/check_mpi_amr_tree_restart_3d.py tools/check_mpi_amr_tree_application_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Build and MPI ABI 40 executable audits passed. CTest **2/2 passed**, 63.02 s
total: selected static application 6.87 s, restart 56.12 s; no failures/skips.
Python syntax, provenance, project contract and whitespace checks passed.
Existing compiler warnings remain.

The restart fixture uses the branched three-level elementary hotspot with
`cfl=0.002`: five steps to 1e-9 s, checkpoint after step 1 at
2.1753230644736716e-10 s. For both fixed and selected mechanisms:

- Two-rank uninterrupted output equals two-rank checkpoint/three-rank restart
  CSV and budget bytes.
- MPI checkpoint/serial continuation and serial checkpoint/one-rank continuation
  produce the same final output. Serial and MPI checkpoint bytes also match.
- Checkpoint stop preserves existing final outputs.
- Changed physics, truncated footer, finite but inconsistent boundary ledger,
  checkpoint/output alias and malformed namelist are rejected, preserving outputs.

The static application regression additionally covers root-only topology with
empty ranks, and now rejects the still-unsupported adaptation group. No full
configuration matrix, installed-prefix audit or external/physical validation
was run. The restart test exercises the selected elementary fixture, not all
supported bundles or chemistry integrators.

Logs: `build/amr-tree-mpi-restart-checks/{build,ctest}.log`.
Per-run inputs/checkpoints/logs/outputs:
`build/selected-tree-mpi-debug/tests/mpi_amr_tree_restart_3d/{fixed,selected}/`.

## Limits

Only static synchronous topology is connected. MPI tree regridding and subcycling
remain open. Root checkpoint gathers have no scaling qualification. The serial
writer's possible partial file after an I/O failure and compiler/runtime-sensitive
context text remain unchanged. Lexical path checks do not resolve symlinks or
hard links.
