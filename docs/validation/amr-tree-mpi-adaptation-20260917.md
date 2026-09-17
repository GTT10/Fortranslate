# Automatic MPI tree adaptation and restart — 2026-09-17

## Source and scope

Base: `2051e700a37fb8851f218e3fc1a28a43f91a1c1b`, local `main`.
Connects fixed/selected MPI application adaptation namelists, distributed tagging,
sparse remapping, adaptation context/counters and changed-rank geometry restoration.
Numerical kernels, mechanism sources, thresholds/tolerances and golden hashes are
unchanged. Scheduled and automatic modes remain mutually exclusive.

Final SHA-256 (driver indentation cleaned after compilation):

- `src/driver/mpi_amr_patch_tree_application_3d_mod.F90`:
  `d1906169fcbfa48fc22af7e58103844697e8439f19417805bb6bd47185803932`
- `tools/check_mpi_amr_tree_adaptation_3d.py`:
  `721c05940dfb6fae0f0e25979999d4259995ff320e7a19633978daa5ff514a3d`
- `tools/check_mpi_amr_tree_application_3d.py`:
  `1315e55825cad029a136e84810544d41de3fd925dbef9e100e2022b8a3d65703`

## Completed checks

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, existing clean system-MPI configuration
`build/selected-tree-mpi-debug`; compiler/backend settings unchanged.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --parallel 4 --target pelef_mpi_amr_tree_reactive_3d test_selected_elementary_mpi_amr_tree_3d test_selected_elementary_amr_tree_3d pelef_amr_tree_reactive_3d
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_mpi_amr_tree_(adaptation_3d|schedule_3d|selected_application_3d)$' --output-on-failure
python3 -m py_compile tools/check_mpi_amr_tree_adaptation_3d.py tools/check_mpi_amr_tree_application_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Build and executable MPI ABI 40 audits passed. CTest **3/3 passed**, 86.18 s total:
selected static application 6.76 s, scheduled regrid/restart 51.07 s, automatic
adaptation 28.32 s; no failures/skips. Python syntax, provenance, project contract
and whitespace checks passed. Existing compiler warnings remain.

The automatic fixture starts from the root-only 4-cubed geometry with an off-center
hotspot, chemistry/transport disabled, CFL 0.001, thresholds 1150/1150 K, ratios
2/2, unit tiles, zero buffer and a 64-node budget. It generates 10 nodes on three
levels and performs seven steps/adaptation checks to 1e-9 s. Both mechanism paths
report zero normalized remap error for this fixture.

For fixed and selected elementary mechanisms, the checker verifies exact serial
versus 2-rank CSV/budget bytes, 2-rank stop/3-rank continuation, no repeated step-1
adaptation, and unchanged output on checkpoint stop. MPI checkpoint/serial restart
also matches; serial and MPI saved checkpoint bytes are identical. The final
fixed/selected outputs agree. Output exposed-cell volume equals domain volume
within 1e-20 absolute tolerance.

Changed threshold/cadence, a corrupted adaptation-check counter and insufficient
node budget are rejected while preserving output. Both counter-rejection logs
explicitly report `MPI tree restart adaptation/schedule position mismatch`.
The static regression now checks rejection of the still-unsupported subcycling
namelist. Scheduled regression retains reactive/transport evolution and event
continuation coverage; the automatic test itself is hydro-only.

Logs: `build/amr-tree-mpi-auto-checks/{build,ctest}.log`.
Automatic run artifacts:
`build/selected-tree-mpi-debug/tests/mpi_amr_tree_adaptation_3d/{fixed,selected}/`.
Other artifacts remain in the `mpi_amr_tree_schedule_3d/` and
`mpi_amr_tree_selected_application_3d/` sibling directories.

## Remaining scope

MPI coarse-time donors, level substeps and recursive subcycling are not connected.
Hot-cell metadata replication, quadratic tagging/remap planning and root
checkpoint/output gathers have no production scaling qualification. No full
configuration matrix, installed-prefix audit, external solver parity or physical
validation was run. Existing lexical-path and checkpoint-I/O limits remain.
