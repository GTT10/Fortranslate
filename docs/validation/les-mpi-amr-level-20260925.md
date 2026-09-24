# MPI level-local LES integration — 2026-09-25

Implementation `5243741`, based on `37b91ea`.
`advance_mpi_amr_tree_les_level_3d` advances one depth with owner-local SSPRK2
and returns time-integrated face fluxes only for locally owned active nodes.
It uses the shared distributed temporal donor sampler, including wide Dynamic
halos, and the existing static/Dynamic wall and pressure kernels.

Model/filter, numerical controls, endpoint presence, depth and physical inputs
are agreed collectively. Temporal fractions are checked by the halo exchange
for each stage. State, temperature, diffusion maximum and integrated fluxes
publish only after both stages and final recovery succeed. Nonowners remain
unallocated; inactive owned levels remain unchanged. No whole-tree state gather
is introduced. Coarse endpoint fields pack conserved state and a temperature
guess, matching the serial contract.

The final temperature guess is the first-stage recovered temperature, matching
serial LES. The MPI maximum diffusion bound includes all active owners, with
zero local contribution on empty ranks. Recursive scheduling and coupled
application subcycle dispatch remain separate work; this level API does not
perform interval reflux by itself.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug. Fixed tests-enabled runtime/test
and selected H2O2 tests-disabled runtime built. Selected execution was not run.
CTest `spray_les_mpi_amr_level_2`: **1/1 passed**, 21.48 seconds, total 21.52.

Smagorinsky, WALE, Vreman and Dynamic exercise depths 0/1/2, refinement ratios
2 and 3, changing coarse endpoints, PCM/MC donors, periodic and moving
isothermal walls, and reservoir inflow. Concentrated ownership includes an
empty rank; distributed ownership puts coarse and fine donors on different
ranks. State, temperature, integrated x/y/z fluxes and maximum diffusivity
match serial exactly. Inactive state and temperature are unchanged.

Excessive dt rejects without changing gas, temperature, flux allocation/values
or the diffusion diagnostic. A second-stage fraction disagreement also rejects
and preserves all outputs after first-stage work has occurred. Pressure-target
and filter-disagreement cases were not separately added to this level test.

No failures, skips or tolerance changes occurred. Mechanism provenance,
Python syntax of the contract checker, project contract and whitespace passed.
This is distributed operator evidence, not recursive MPI integration, external
solver parity, physical validation or a full configuration matrix.

## Commands and artifacts

Logs: `build/les-mpi-amr-level-checks/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_mpi_amr_level --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_mpi_amr_level_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
