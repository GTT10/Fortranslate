# Stage-time gas halos connected to EB ghost evolution — 2026-09-27

Implementation `53441199`, based on `cf46aaed`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`advance_mpi_eb_ghost_stage` fills complete conserved/temperature halos with
the existing MPI EB gas halo operator, then calls the integrated ghost stage.
Ancestor start/end states, temperatures and time fractions are forwarded.
`restrict_snapshot=.false.` preserves the caller's stage states. Target volume
fractions and centroids are caller-supplied EB metrics. The gas states are
read-only; halo-fill failure returns before any particle/source/ledger update.

This connects stage-time gas data with kick, drift, reaction transfer and
redistribution. Recursive clocks, physical/virtual/ghost accounting, breakup
and application dispatch are still unconnected. The underlying EB ghost
redistribution retains its nonperiodic restriction.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. `spray_mpi_eb_sample_2` and `_4` passed
(2/2, 0.60 s and 0.59 s; total 1.23 s).

The existing two-level planar EB fixture now evolves one water ghost in the
fine patch's first outer x cell. Root gas endpoints are U and 2U at 900 K.
Fractions 0, 0.25 and 1 select the expected parent gas states. A 1e-8 s kick
with zero drift matches the frozen-gas reference using those states: velocity
within 1e-12, temperature within 1e-10 and mass within relative 1e-14. The
fluid-volume integral of the routed source matches the extensive reaction
within relative 1e-12 (denominator floor 1e-30). Fraction 1.1 rejects while
preserving components, source fields and sentinel wall/boundary ledgers.
Existing EB interpolation and rollback assertions also execute.

Provenance preflight, target build, Python syntax, project contract and
whitespace checks passed. No build/test failures or skips occurred. The new
wrapper is tested without drift here; this is not a temporal-order check,
recursive application run, full matrix, installed audit, external solver parity
or physical validation.

## Reproduction

Artifacts: `build/spray-ghost-temporal-checks/`, including copied CTest details.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_sample --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sample_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
