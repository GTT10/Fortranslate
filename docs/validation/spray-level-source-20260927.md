# Level-bounded particle source deposition — 2026-09-27

Implementation `3f9026a7`, based on `7c0f3397`, local `main`, version 0.246.0.
PeleC reference `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance preflight SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Source semantics

`deposit_mpi_amr_tree_spray_field` accepts optional `maximum_level` for
constructing temporary level-view source fields. The capped regular/EB cloud
uses only cells exposed at or below the requested ceiling. Contributions are
normalized by cell/fluid volume and accumulated with the existing ordered
MPI exchange and overflow guards. Exported `gas_change` retains independent
increments rather than subtracting the background field.

An explicit ceiling suppresses parent restriction, including when it equals
the finest level. Otherwise deeper unchanged states would erase source deposits
in covered coarse cells. Cells above the ceiling remain unchanged. Covered
ancestors are not refreshed: the caller must treat this as raw source-field
construction. Omission retains the existing physical composite deposition and
restriction. Presence and effective value of the ceiling agree across ranks
before any update. The gas/EOS wrapper does not expose this option.

This is a prerequisite for separate real/virtual/ghost level source fields.
It does not transfer/synchronize those representations or prevent double-counting
by an incorrectly assembled caller. Fine ghost interpolation, source transfer
and recursive spray application integration remain open. No full solver parity
or physical validation follows from this field-level test.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6 through `/usr/bin/mpifort` and `/usr/bin/mpiexec`.

Six CTests passed, 3.03 s total:
`spray_mpi_amr_deposition_2`, `_4`, `spray_mpi_eb_sample_2`, `_4`,
`spray_mpi_eb_periodic_sample_2`, `_4`.

The regular three-level fixture deposits at ceilings zero and one and checks
exported increments reconstruct the updated field within 1e-12, deeper fields
remain exactly unchanged, and the capped-view volume integral equals summed
input contributions within 1e-12. A zero receive capacity rejects collectively
and preserves the field. Existing default composite serial/MPI and rollback
checks continue to execute.

The EB fixture adds coarse source deposition to its existing sampling checks.
Both periodic and nonperiodic cases integrate the exported source with coarse
fluid volumes and recover the prescribed contribution of 2 within 1e-13.
Fine fields remain exactly unchanged and fine exported increments are zero.

Provenance, target builds, Python syntax, project contract and whitespace checks
passed. No build/test failures occurred. Existing compiler warnings remain.
No full matrix, selected-mechanism build, installed audit, recursive application
run or physical validation was performed.

## Reproduction

Artifacts: `build/spray-level-source-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_deposition test_spray_mpi_eb_sample --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(amr_deposition|eb_(periodic_)?sample)_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
