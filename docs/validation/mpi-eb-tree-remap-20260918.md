# Sparse MPI EB remapping and scheduled regrids — 2026-09-18

Implementation: `0beaf02fcf95e843ec68f4f5416996dc57647944`, based on
`b82451b`. Module `mpi_eb_tree_remap_3d_mod` is part of `pelef_mpi_support`;
the shared application now dispatches prescribed MPI topology changes.

## Implemented behavior

`remap_mpi_eb_tree_3d` accepts old/new topologies and owner maps, a finest
geometry covering their deeper level, and owner-local U/T. It first averages a
private old snapshot. New patch owners request each destination's deepest
available old donor at the same or a coarser level. Averaged parents supply
coarsening; exact same-level copies or ancestor PCM supply overlap/refinement.
Final averaging precedes collective publication of geometry, U and T together.
Caller topology/owner metadata changes only after success.

The API requires a common root and unchanged refinement ratios on shared
levels. Coarsened finest geometry must match all old patch metrics; new metrics
must agree across ranks. No full state gather is used. Sparse requests include
inactive cells to preserve their payloads; as required by the reused MPI cell
exchange, fetched U/T must be finite even there. EOS validation applies to
fluid destinations. This is not moving-geometry or arbitrary-ratio remapping.

The application now accepts prescribed regridding, lands exactly on event
times, chooses new whole-patch owners and dispatches distributed remapping.
Its geometry-only tree metadata are rebuilt after success. Existing scheduled
checkpoint envelopes support changed-rank and serial interchange before/after
events. Saving/CSV remains root-gathered. Temperature-driven adaptation and
subcycling remain unsupported by the MPI entry point.

## Completed focused verification

Fresh GNU Fortran 13.3.0 Debug/Ninja build; compiler `/usr/bin/gfortran`, MPI
wrapper `/usr/bin/mpifort`, launcher `/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake -S . -B build/mpi-eb-tree-remap-debug -G Ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_Fortran_COMPILER=/usr/bin/gfortran -DMPI_Fortran_COMPILER=/usr/bin/mpifort -DMPIEXEC_EXECUTABLE=/usr/bin/mpiexec -DPELEF_ENABLE_MPI=ON -DPELEF_ENABLE_TESTS=ON
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/mpi-eb-tree-remap-debug --target test_mpi_eb_tree_remap_3d pelef_reactive_eb_tree_3d pelef_mpi_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/mpi-eb-tree-remap-debug -R '^(unit_mpi_eb_tree_remap_3d_[23]|regression_mpi_eb_tree_scheduled_3d)$' --output-on-failure --output-log build/mpi-eb-tree-remap-checks/ctest.log
python3 -m py_compile tools/check_mpi_eb_tree_application_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**3/3 tests passed**, no failures/skips, total **26.32 s**. Two/three-rank
unit tests took 1.20/1.21 s; the scheduled application test took 23.89 s.

Unit fixtures cover three plane orientations and four/five-node trees, with
patch movement, full coarsening and refinement back to three levels. Ownership
changes at each stage, including a root-only stage with empty ranks. Remapped
U/T match serial exactly. The **36** MPI inventory comparisons have maximum
normalized residual **7.12457338238146e-15**, below the unchanged 3e-14 bound.
Invalid source EOS and changed geometry preserve geometry/U/T.

The fixed explicit application changes **5 nodes -> 1 -> 4** at 1.5e-9 and
2.5e-9 s. Serial/two/three-rank uninterrupted runs and serial/MPI restart
directions produce identical final state CSVs. Two-to-three-rank continuation
before and after both events also matches. Same-rank restart matches the budget
CSV exactly; changed-rank budgets agree within 3e-14 scaled by initial inventory.
Independent volume integrals and ledger checks pass with maximum normalized
residual **7.05347627465518679e-15**. Existing corrupt-record, changed cap,
automatic-mode/CFL policy and unsupported-subcycling rejection checks protect
sentinel outputs.

Configure/build, provenance, Python syntax, project contract and whitespace
checks pass. Compiler warnings remain. This does not qualify automatic-plus-
scheduled integration independently, selected MPI, temperature adaptation,
subcycling, installed binaries, scaling, a full matrix or physical validation.

Logs: `build/mpi-eb-tree-remap-checks/{configure,build,ctest}.log`.
Runtime artifacts: `build/mpi-eb-tree-remap-debug/tests/mpi_eb_tree_scheduled_3d/`.
