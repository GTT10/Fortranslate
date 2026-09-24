# AMR reacting transport/hydro Mach inlet composition — 2026-09-25

Implementation `25be739`, based on `5824244`.

## Change

Serial and MPI `R/2 - T/2 - H - T/2 - R/2` drivers forward optional Mach and
inlet-pressure arrays through common or recursive hydro and both molecular
transport half steps. Inlet temperature uses the existing argument. Chemistry
remains unchanged; candidate state and diagnostics are published only after the
full composition succeeds. Spray coupling and application input/timestep
integration remain pending.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug. Fixed MPI
test and selected H2O2 spray-tree runtime dependencies built; selected execution
was not run. `unit_mpi_amr_patch_tree_split_3d_3` passed: **1/1 tests**, no failures
or skips, 60.91 s (CTest total 60.95 s).

New Mach 0.2 cases enable chemistry and all molecular terms. The explicit/common
case extrapolates pressure; the implicit/recursive case prescribes 210000 Pa and
uses moving isothermal walls. State, temperature, boundary/reflux diagnostics
and minimum limiter match serial exactly. The third rank owns no nodes in the
Mach cases. An invalid hydro solver after preceding reaction/transport work
preserves the original state and diagnostics; inconsistent operator flags reject
collectively. Existing reservoir, periodic, wall and outlet cases remain.
Provenance, Python syntax, project contract and whitespace passed. This is
focused coupled parity, not physical ignition or nozzle validation.

## Commands and artifacts

Logs: `build/mach-inlet-amr-split-checks/`, containing `build-fixed.log`,
`build-selected.log` and copied `ctest.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_split_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_amr_patch_tree_split_3d_3$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Environment overrides are host-local workarounds, not application defaults.
