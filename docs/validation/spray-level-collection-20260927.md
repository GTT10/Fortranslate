# Collection of synchronized physical level populations — 2026-09-27

Implementation `2db18429`, based on `853814b5`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`spray_level_population` groups particle, component and node arrays.
`collect_mpi_amr_level_particles` validates physical level snapshots, removes
dead entries, exchanges live entries by ID modulo MPI size, and uses the
existing sorted exchange to reject duplicate live IDs globally. It then
migrates the private population to current finest geometric owners, including
optional EB validation. Stale node indices are ignored. Caller outputs publish
together only after all operations succeed; the input snapshots remain read-only.

Local live input and temporary ID-owner bins retain the existing one-million
particle bound. The caller supplies only physical populations at synchronized
times; virtual and ghost copies must remain separate. Clock management and
integration of split/advance/collect into the recursive application remain open.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. Four tests passed in 1.98 s total:

- `spray_mpi_eb_tree_2`: 0.48 s.
- `spray_mpi_amr_tree_2`: 0.48 s.
- `spray_mpi_eb_tree_4`: 0.49 s.
- `spray_mpi_amr_tree_4`: 0.49 s.

The regular three-level fixture moves the selected physical populations across
coarse/fine boundaries before collection. It verifies current owners/positions,
unchanged components and exactly one occurrence of every physical ID globally.
Capacity zero preserves saved output IDs/components. A duplicate ID is then
inserted at a different position; the ID exchange rejects it and preserves the
previous outputs. The four-rank case places those duplicate positions on
different geometric owners. Existing virtual, ghost and regular/EB migration
assertions also execute.

Provenance, build, Python syntax, project contract and whitespace checks passed.
No build/test failures or skips occurred. Dead-entry filtering and the new
collection API's EB path were compiled but not independently exercised here.
No full matrix, installed audit, recursive application, external solver parity
or physical validation was performed.

## Reproduction

Artifacts: `build/spray-level-collection-checks/`, including CTest details.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_tree test_spray_mpi_eb_tree --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(amr|eb)_tree_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
