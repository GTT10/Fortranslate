# Explicit sparse MPI hierarchy repartition — 2026-09-17

Source: parent `891d667` plus this change. Logs and source hashes:
`build/amr-mpi-repartition-checks/`.
Reused GNU Fortran 13.3.0 Debug / Open MPI 4.1.6 build with bounds checks.

The sparse distribution initializer now optionally accepts nonuniform coarse
slab counts. Each rank retains at least one coarse column; fine ownership is
the intersection with the patch's parent cells. The explicit repartition API
collectively validates the old layout and requested counts, constructs a
private target layout and transfers coarse/fine states plus temperatures with
`MPI_Alltoallv`. Payload counts are checked against the default MPI integer
range before allocation or communication. Empty peer transfers are skipped.
There is no root field aggregation and no arithmetic on the transferred values.
Hierarchy and ownership commit together after both levels pass validation.

The conservative remap now preserves the current coarse counts when creating
new fine ownership. Its parent-cell locality assumption therefore remains
valid after explicit repartitioning. Automatic work estimates and application
load-balancing policy are not yet connected; application defaults are unchanged.

The existing four/eight-rank remap test now interleaves explicit repartitioning
with movement, growth, shrinkage and full coverage at refinement ratios two
and three. Four ranks alternate heavily left/right weighted ownership; eight
ranks retain the minimum one coarse column per rank. Global states and
temperatures remain exactly equal before/after transfer and after serial/MPI
remap. The remap preserves requested coarse counts. Invalid and one-rank
mismatched valid counts reject collectively without changing either level or
coarse/fine ownership. Existing geometry and inlet rollback checks remain.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target test_mpi_amr_regrid_3d --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^mpi_amr_regrid_3d_np(4|8)$' --output-on-failure
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 2/2 tests passed, no failures or skips (four ranks 1.40 s, eight
ranks 1.08 s; total 2.51 s). Provenance, project contract and diff checks
passed. No Python source changed. The build succeeded with real-equality
warnings. Only the final run after the ownership assertions and empty-peer
skip was added is counted here.

This is focused data-transfer/remap evidence. Time stepping and checkpoint
application behavior with automatically chosen uneven partitions, production
scaling, external solver parity and physical validation were not qualified.
General patches/levels and distributed I/O remain open.
