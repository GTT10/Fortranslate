# MPI AMR particle physical boundaries — 2026-09-18

Implementation: `07e7879142908185cb44d79deeb80abd48dbae6f`, based on
`defbc3a`. `advance_mpi_amr_mixture_coupling` now accepts face types and an
optional integrated boundary ledger. It reuses existing mixture-property and
parcel-boundary kernels: periodic wrapping, removal at inflow/outflow, and
specular reflection at gas wall face types. This is not wall impact/film or
liquid thermal-wall modeling.

Face types agree collectively before stencil queries. Per-trial boundary
inventory is private, accepted substeps accumulate it, and a global reduction
publishes only with successful final migration. Failed intervals preserve the
incoming ledger along with gas, particle/component arrays and IDs. Gas feedback
uses pre-step positions; exited parcels have zero component inventory and are
removed from final survivors. Existing periodic calls remain supported.

## Completed focused verification

Reused GNU Fortran 13.3.0 Debug/Ninja fixed tests-enabled and selected H2/O2
tests-disabled builds, `/usr/bin/gfortran`, Open MPI `/usr/bin/mpifort` and
`/usr/bin/mpiexec`. Compiler/MPI/backend configuration did not change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_boundaries test_spray_mpi_amr_step --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_(boundaries_(2|4)|step_2)$' --output-on-failure --output-log build/spray-amr-boundary-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**3/3 passed**, no failures/skips: two-rank boundaries **0.53 s**, existing
two-rank periodic interval **0.74 s**, four-rank boundaries **0.54 s**,
**1.84 s** total. The three-level fixture removes one parcel at the high-x
outflow, reflects another at the low-x wall and retains an interior parcel.
Independent gas-plus-liquid integrals agree with the returned mass, momentum,
energy and species ledger within 3e-12 of inventory/ledger scale. Explicit
assertions check outgoing represented mass 1e-10 kg and net x momentum change
1e-7 kg m/s, final count and reflected velocity/direction.

Substep exhaustion after a boundary event, rank-disagreed face types and late
receive-capacity failure all preserve original gas, temperatures, parcels,
component masses, IDs and a sentinel boundary ledger. The boundary fixture
disables drag/heat/evaporation to isolate boundary bookkeeping; the retained
periodic regression exercises those models and KH children separately.

The selected library compiled successfully; this is compile evidence, not a
selected physical-boundary executable run. All build/test attempts passed.
Provenance, Python syntax, project contract and whitespace checks passed.
Compiler warnings remain. No Python source changed. Reacting-step/application
physical-boundary integration, checkpoint boundary ledger, wall films, EB
collisions, external parity and physical validation remain outside this check.

Logs: `build/spray-amr-boundary-checks/{build,build-test,build-selected,ctest}.log`.
Executables: `build/spray-amr-reactive-debug/tests/{test_spray_mpi_amr_boundaries,test_spray_mpi_amr_step}`.
