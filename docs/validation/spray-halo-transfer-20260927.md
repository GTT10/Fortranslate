# Reaction transfer inside the particle halo transaction — 2026-09-27

Implementation `69f3a8e8`, based on `313fdf89`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`kick_mpi_amr_mixture_halo` accepts the optional group `target_level`,
`level_source`, `source_remainder`. Each queried node must have the selected
level. After private kick, optional drift and deposition, the existing source
transfer sums the reaction into a fresh owner-local valid-cell field and
returns unmatched contributions separately. Optional `route_to_coarse` enables
conservative ancestor routing; `maximum_received` bounds exchange capacity.
EB transfer requires both destination `eb_geometry` and source `fluid_fraction`.

All source outputs, particles, components and optional drift ledgers publish
only after successful transfer. The original raw source halo remains an output.
This connects previously separate particle and source-transfer operations;
it does not apply the source to gas, schedule recursive levels, migrate moved
populations, or deduplicate physical/virtual/ghost inventories.

## Completed verification

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. `spray_mpi_amr_step_2` and `_4` passed
(2/2, 0.82 s and 0.78 s; total 1.65 s).

The fine-halo kick fixture now requests integrated transfer at level one with
coarse routing. Its particle lies just beyond the fine patch's valid x extent.
The coarse owner receives a nonzero source. The volume integral over coarse
and fine valid fields matches the extensive kernel reaction within relative
1e-12 (1e-30 denominator floor), and all remainder entries are zero.
Repeating with receive capacity zero fails after private particle updates;
components, checked velocity/temperature and raw/routed/remainder source fields
retain their saved values exactly. Existing kick, periodic drift, outflow,
composite coupling and capped-view assertions execute in the same tests.

Provenance, target build, Python syntax, project contract and whitespace checks
passed. No build/test failures or skips occurred. Integrated EB transfer and
combined drift-plus-transfer were not separately exercised here. No full matrix,
installed audit, recursive application, external solver parity or physical
validation was performed.

## Reproduction

Artifacts: `build/spray-halo-transfer-checks/`, including copied CTest details.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_step --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_step_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
