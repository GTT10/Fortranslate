# Distributed dynamic LES time integration — 2026-09-18

Implementation `6fa53df41a264e79bc44b99e121f9e7492e4f32e`, based on `d2b6b05`.
`advance_mpi_dynamic_les_ssprk2_3d` advances owner-local periodic x slabs with
stage-wise dynamic coefficient/flux evaluation. All ranks must supply the same
finite dt. Allocation, RHS and EOS acceptance are collective; state/temperature
publish only after the second recovery succeeds. Each stage uses the global
maximum frozen-coefficient diffusion estimate and the same 0.2 diffusion-number
gate as serial. Zero dt still validates the distributed model/state, then returns
without mutation. The heuristic timestep limitation remains; this is not a
nonlinear stability proof. MPI spray/application dispatch and nonperiodic/AMR
paths remain open.

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; unchanged configurations.
Fixed tests-enabled `build/spray-amr-reactive-debug` built `test_les_mpi_dynamic`.
Selected H2O2 tests-disabled `build/spray-amr-reactive-selected-debug` built
`pelef_selected_spray_runtime`.

Focused CTest **2/2 passed**, 1.04 seconds, no failures/skips. Two/four ranks use
a 5-by-3-by-2 variable-density/temperature/velocity field with anisotropic spacing,
unequal slabs and one-cell local extents. Box and Gaussian configurations produce
state and temperature exactly equal to the serial SSPRK2 result. Global component
conservation passes relative `1e-12` with unit floor. Zero dt is exactly unchanged;
a diffusion-number 0.5 step and a last-rank dt mismatch reject with exact state
and temperature preservation. Existing distributed coefficient/flux parity and
configuration/EOS rejection checks remain. Provenance, Python syntax for the
contract checker, project contract and whitespace checks passed.
Selected runtime was compiled, not exercised. No full matrix, external solver
parity, large-scale MPI or physical validation ran.

## Commands and artifacts

Logs: `build/les-mpi-dynamic-step-checks/{provenance,build,build-selected,ctest}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_dynamic --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_mpi_dynamic_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
