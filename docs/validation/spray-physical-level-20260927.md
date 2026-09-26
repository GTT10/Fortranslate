# Disjoint physical particle populations by AMR level — 2026-09-27

Implementation `a7bcb490`, based on `8f8ee790`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`create_mpi_amr_level_particles` selects particles whose finest containing
node has exactly the requested level. It shares validated selection/exchange
code with virtual-particle creation, whose strictly-finer selection remains
unchanged. MPI agrees on selection kind and level before processing.
The new routine returns owner-local parcels, components and node indices
matching the exchanged ordering. All outputs remain private until success.
Optional nodal EB geometry follows the existing migration/solid checks.

The source is read-only. Level outputs are snapshots partitioning physical
particles, not additional physical inventory. Recursive clocks, advancement
and final merging remain to be connected; this is not a completed population
scheduler or spray application.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. Four tests passed in 2.01 s total:

- `spray_mpi_eb_tree_2`: 0.48 s.
- `spray_mpi_amr_tree_2`: 0.48 s.
- `spray_mpi_eb_tree_4`: 0.50 s.
- `spray_mpi_amr_tree_4`: 0.49 s.

The regular three-level fixture selects each level in turn. A global per-ID
counter proves every input physical ID appears exactly once across the outputs.
It checks exact level/owner, position, velocity, temperature and component
inventory. The fully covered root has no physical population. Capacity zero
rejects and preserves saved IDs/count, components and nodes; the input IDs and
components remain unchanged. Existing virtual, ghost and regular/EB migration
assertions also execute.

The new exact-level API's EB path was compiled but not independently exercised.
Provenance, build, Python syntax, project contract and whitespace checks passed.
No build/test failures or skips occurred. No full matrix, installed audit,
recursive application, external solver parity or physical validation was run.

## Reproduction

Artifacts: `build/spray-physical-level-checks/`, including CTest details.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_tree test_spray_mpi_eb_tree --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(amr|eb)_tree_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
