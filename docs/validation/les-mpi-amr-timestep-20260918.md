# MPI AMR LES timestep — 2026-09-18

Implementation `599a98e9ab13400df34941471d0a53604b7da007`, based on `c58ba4c`.

`compute_mpi_amr_tree_les_timestep_3d` selects a common diffusion interval from
owner-local fields and distributed ghosts. Private covered-cell averaging and
EOS recovery precede evaluation, matching the initial integration stage. Every
patch uses its own spacing and the shared static LES RHS; MPI reduces interval
by MIN and diffusivity by MAX. Empty owners contribute the interval cap and zero
diffusivity. No global state gather is introduced. CFL is finite in `(0,0.2]`,
with the same `1e6 s` cap as the serial selector. Model, coefficients, CI, CFL,
spacing and boundary data agree collectively. Outputs publish only after all
owners and reductions succeed; caller state is read-only.

This is an initial frozen-gradient estimate. Both integration stages still
check their own interval and EOS constraints; application retry remains needed.
Subcycling, dynamic LES, EB and application dispatch are not part of this change.

## Completed verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; existing configurations unchanged.
Fixed tests-enabled `test_les_mpi_amr_tree` and selected H2O2 tests-disabled
`pelef_selected_spray_tree_runtime` built. Selected runtime was not run.
CTest **2/2 passed**, 27.36 seconds, zero failures/skips:
`spray_les_mpi_amr_tree_2` (14.14 s), `spray_les_mpi_amr_tree_4` (13.18 s).

For all three closures with CI 0.1, root/six-node trees, concentrated/distributed
owners and periodic/outflow/moving-isothermal plus pressure boundaries, MPI dt
and maximum diffusivity match the serial selector within `2e-12` relative.
Rank-dependent CFL values reject collectively while preserving dt/diffusivity
sentinels. Existing owner-state/temperature/ledger parity and rollback checks
also pass. Step size remains `1e-9 s`, so this is not acceptance at the selected
upper bound or a scaling qualification. No tolerances changed.
Provenance, Python syntax, project contract and whitespace checks passed.
No failed builds/tests, full matrix, external solver parity or physical
validation occurred.

## Commands and artifacts

Logs: `build/les-mpi-amr-timestep-checks/{provenance,build,build-selected,ctest}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_amr_tree --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_mpi_amr_tree_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
