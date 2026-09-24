# MPI AMR Mach inlet ghost exchange — 2026-09-25

Implementation `622ce86`, based on `2098aad`.

## Change

The collective AMR ghost entry forwards optional species, Mach, inlet temperature
and pressure to every sampler pass: option validation, donor planning and final
state construction. After collective validation, all ranks must agree on the
effective Mach/pressure arrays and active inlet temperatures before exchanging
donors. Empty ranks participate in validation and agreement. Output publication
remains transactional after collective success. Existing donor routing, temporal
endpoints, PLM interpolation and depth selection are retained.

AMR gas/LES stepping and application input are still pending; this change connects
the distributed ghost entry only.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug.
Fixed MPI unit target and selected H2O2 spray-tree runtime built. Selected
execution was not run. **2/2 CTests passed**, no failures or skips, total 37.42 s:

- `unit_mpi_amr_patch_tree_ghost_3d_2`: 18.67 s.
- `unit_mpi_amr_patch_tree_ghost_3d_3`: 18.71 s.

New checks exercise all six inlet faces, seven ghost layers, extrapolated and
prescribed pressure, MC interpolation and two coarse temporal fractions in a
three-level tree. Owned output agrees exactly with the serial sampler. The third
rank owns no nodes. Mach, temperature and pressure disagreement, and invalid
temperature on the last rank reject collectively and preserve all previous
outputs. Existing PCM/PLM, temporal, physical-wall and rollback checks also pass.
Mechanism provenance, Python syntax, project contract and whitespace passed.
This is focused ghost-exchange evidence, not physical nozzle validation.

## Commands and artifacts

Logs: `build/mach-inlet-mpi-amr-ghost-checks/`, including `build-fixed.log`,
`build-selected.log` and copied `ctest.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_ghost_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R 'mpi_amr_patch_tree_ghost_3d' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Environment overrides are host-local workarounds, not application defaults.
