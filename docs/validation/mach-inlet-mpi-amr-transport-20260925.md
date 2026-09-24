# MPI AMR Mach inlet molecular transport — 2026-09-25

Implementation `fb92c88`, based on `412277e`.

## Change

Common-step and per-level MPI molecular transport forward optional Mach,
temperature and inlet pressure through both SSPRK2 ghost fills. Active conditions
are validated and agreed collectively before reservoir handling. The recursive
driver forwards them through every interval. Existing viscous wall handling,
conservative limiting, restriction, reflux and transactional output remain.
AMR LES, coupled dispatch and application input/timestep integration remain pending.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug. Fixed
MPI unit targets and selected H2O2 spray-tree runtime dependencies built.
Selected execution was not run. **2/2 CTests passed**, no failures or skips,
total 115.97 s:

- `unit_mpi_amr_patch_tree_transport_3d_3`: 9.31 s.
- `unit_mpi_amr_patch_tree_transport_subcycle_3d_3`: 106.63 s.

New Mach 0.2 cases use extrapolated or 120000 Pa inlet pressure, all molecular
terms and existing physical walls. State, temperature, boundary/reflux diagnostics
and limiter values match serial exactly; recursive counts also match. The third
rank owns no nodes in the Mach cases. Inlet-pressure disagreement on the last
rank preserves state and diagnostics. Existing transport-property failure,
timestep disagreement and outlet-pressure guards remain included. Mechanism
provenance, Python syntax, project contract and whitespace passed. No tolerance
or timeout changed. This is focused numerical parity, not physical validation.

## Commands and artifacts

Logs: `build/mach-inlet-mpi-amr-transport-checks/`, containing `build-fixed.log`,
`build-selected.log` and copied `ctest.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_transport_3d test_mpi_amr_patch_tree_transport_subcycle_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_amr_patch_tree_transport(_subcycle)?_3d_3$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Environment overrides are host-local workarounds, not application defaults.
