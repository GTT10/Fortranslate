# Caller-assembled predictor halos for MPI EB hydro — 2026-09-27

Implementation `16c20b57`, based on `3c25ed34`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`compute_sdc_mpi_eb_tree_hydro_3d` accepts optional `predictor_halo`. The normal
predictor fill still validates source fields and temporal metadata. If supplied,
the complete three-cell halo must match owner allocation and shape, contain
finite values, and exactly match the normal predictor at all valid fluid cells.
Spatial lower-bound offsets permit packed or physical layouts. MPI agrees on
argument presence. After collective validation, the supplied halo replaces the
normal halo; values are not added a second time.

This permits ghost-region particle sources to reach CTU prediction without
being reconstructed from valid-cell sources. The caller must assemble the
complete predictor, including all relevant source families. Wiring such an
assembly into the recursive particle callback remains open.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. `unit_mpi_eb_tree_ctu_hydro_3d_3` passed
(1/1, 22.81 s; CTest total 22.85 s).

The existing three-level EB test runs two ownership layouts: distributed and
all patches on rank zero. Supplying the normally filled predictor halo gives
exactly the baseline hydro rates. Changing only its lower-x ghost momentum
source changes the hydro flux divergence. Changing valid-cell values then
rejects and preserves all checked rate and x/y/z flux arrays. Existing serial
parity, selected-level, temporal-source and failure checks also execute.
The new fixture uses packed halo bounds; physical lower-bound support is
implemented but was not separately exercised here.

The initial build failed because a new test MPI call used an undeclared status
variable. A local declaration fixed it; the initial failure log is retained.
Final build, provenance, Python syntax, project contract and whitespace checks
passed. No test failures or skips occurred. No full matrix, installed audit,
recursive spray application, external solver numerical parity or physical
validation was performed.

## Reproduction

Artifacts: `build/sdc-predictor-halo-checks/`, including final CTest details.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_ctu_hydro_3d_3$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
