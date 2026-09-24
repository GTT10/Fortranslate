# MPI AMR Mach inlet hydro — 2026-09-25

Implementation `ce4f6f8`, based on `c4f55ad`.

## Change

Distributed common-step and per-level hydro accept optional Mach and inlet
pressure controls. They validate temperature shape and active Mach conditions
collectively, then agree on Mach, active temperature and pressure before choosing
the reservoir path. Each SSPRK2 ghost fill forwards the controls to the shared
MPI sampler. Recursive hydro forwards them through every level interval.
Candidate state, temperature, fluxes and diagnostics retain transactional
publication. Reservoir tables are required only for ordinary inlet faces.

AMR molecular transport, LES, coupled application dispatch and inlet-aware
timestep selection remain pending.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6 (`/usr/bin/mpiexec`), Ninja Debug. Fixed MPI
unit targets and selected H2O2 spray-tree runtime dependencies built successfully.
Selected execution was not run. **4/4 CTests passed**, no failures or skips,
total 55.40 s:

- `unit_mpi_amr_patch_tree_hydro_3d_2`: 6.80 s; `_3`: 6.65 s.
- `unit_mpi_amr_patch_tree_subcycle_3d_2`: 20.90 s; `_3`: 21.00 s.

New cases add Mach 0.2 inlet with PCM/extrapolated pressure and characteristic
PLM/120000 Pa inlet pressure. The third rank has no owned nodes in these cases.
Conserved state, temperature, boundary budget and maximum reflux correction match
serial execution exactly. Recursive step counts also match. Mach disagreement
on the last rank rejects collectively with state and diagnostics unchanged;
existing invalid solver, timestep, outlet-pressure and wall-pressure guards pass.
Existing periodic, reservoir and physical-boundary tests remain included.
Mechanism provenance, Python syntax, project contract and whitespace passed.
This is numerical parity evidence, not physical nozzle validation.

## Commands and artifacts

Logs: `build/mach-inlet-mpi-amr-hydro-checks/`, containing `build-fixed.log`,
`build-selected.log` and copied `ctest.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_hydro_3d test_mpi_amr_patch_tree_subcycle_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_amr_patch_tree_(hydro|subcycle)_3d_[23]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Environment overrides are host-local workarounds, not application defaults.
