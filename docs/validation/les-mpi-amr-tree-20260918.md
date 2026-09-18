# Distributed AMR static LES — 2026-09-18

Implementation `e2819511afa68079ab3a9843a7fea6da7851f787`, based on `aa0e3ed`.

`les_mpi_amr_tree_mod::advance_mpi_amr_tree_les_3d` implements the serial static
LES tree split on owner-local fields. Existing distributed validation, ghost
queries, averaging, reflux and EOS APIs provide the mesh operations. No global
state gather is introduced. Model, coefficients, CI, dt, spacing, boundaries,
pressure outlets, coarse limiter and wall/inflow data agree collectively.
Both stages compute local face fluxes, reject excessive diffusion intervals,
synchronize/reflux and recover temperature. Final EOS recovery uses first-stage
temperature guesses like the serial LES path. State/temperature and boundary/
maximum-correction diagnostics commit only after collective success. Runtime
work timing uses the existing patch work tracker.

Supported models are Smagorinsky/WALE/Vreman, optional CI, common timesteps,
regular patch trees and existing physical boundaries. Empty owners are valid.
An MPI timestep selector, subcycling, dynamic LES, EB and application dispatch
remain open. The existing diffusion heuristic is not a general stability proof.

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; configurations unchanged.
Fixed tests-enabled `test_les_mpi_amr_tree` and selected H2O2 tests-disabled
`pelef_selected_spray_tree_runtime` built. Selected runtime was not run.
CTest **2/2 passed**, 18.76 seconds, zero failures/skips:
`spray_les_mpi_amr_tree_2` (9.76 s), `spray_les_mpi_amr_tree_4` (8.96 s).

Each run covers root-only and six-node nested/adjacent hierarchies, three models
with CI 0.1, periodic/outflow/moving-isothermal-wall plus pressure-outlet cases,
MC coarse interpolation, and all-on-rank-zero or round-robin ownership. Nonowner
payloads remain unallocated. Owner-local states match serial to `2e-12` relative
with unit floor; temperature to `1e-8 K`, boundary ledger to `1e-14` absolute,
maximum correction to `1e-12` relative with unit floor. Serial composite
conservation and rollback assertions remain active. Rank-dependent CI and
excessive dt reject collectively with exact field/diagnostic preservation.
The test uses dt `1e-9 s`; it is not a timestep-limit or scaling qualification.

Provenance, Python syntax, project contract and whitespace checks passed.
No failed builds/tests, full matrix, external solver parity or physical
validation occurred.

## Commands and artifacts

Logs: `build/les-mpi-amr-tree-checks/{provenance,build,build-selected,ctest}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_amr_tree --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_mpi_amr_tree_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
