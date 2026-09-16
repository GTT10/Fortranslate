# Distributed MPI temperature tagging — 2026-09-17

Source: parent `cca2e1f` plus this change. Logs and source hashes:
`build/amr-mpi-distributed-tagging-checks/`.
Reused GNU Fortran 13.3.0 Debug / Open MPI 4.1.6 build with bounds checks.

The sparse MPI tagger scans each rank's coarse temperature slab, reduces three
lower and three upper tag coordinates, and constructs the same buffered patch
on every rank. It retains the serial threshold comparison, empty-tag behavior,
buffer clipping and physical/periodic eligibility rules. Geometry, eligibility,
buffer and threshold contracts are checked collectively before reductions.
Invalid local temperatures or mismatched inputs reject collectively with the
old patch and `changed = .false.`.

The application now calls this tagger before gathering any hierarchy state.
If geometry is unchanged, it proceeds directly to stepping. Only an actual
geometry change gathers coarse/fine state for the existing conservative root
remap and redistribution. This removes full-state communication and root
allocation from unchanged tagging steps. It does not make the remap itself
fully distributed, and is not a measured scaling result.

A small MPI test compares serial and distributed tag results with four/eight
ranks: periodic interiors, mixed physical/periodic faces, all physical faces,
threshold equality, cross-rank bounds, clipped buffers, empty and unchanged
tags. It checks collective rejection for one-rank NaN, threshold/buffer or
eligibility mismatch, and invalid face metadata. The boundary-touching dynamic
application regression retains exact serial CSV comparison through four-to-two
rank checkpoint/restart and the existing context rejection probes.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_mpi_amr_tagging_3d pelef_mpi_amr_reactive_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^(mpi_amr_tagging_3d_np(4|8)|regression_mpi_amr_boundary_touching_dynamic_restart_3d)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 3/3 tests passed, no failures or skips (tagging four ranks 0.41 s,
eight ranks 0.45 s; application 57.97 s; total 58.87 s). Provenance, project
contract and diff checks passed. No Python source changed. The build succeeded
with real-equality warnings. Final application artifacts:
`build/spray-mpi-migration-gnu-debug/tests/mpi_amr_boundary_touching_dynamic_restart_3d/run-fx400v7z/`.

Distributed remap, general patches/levels, external solver parity, production
scaling and physical qualification remain separate work. No full configuration
matrix was run.
