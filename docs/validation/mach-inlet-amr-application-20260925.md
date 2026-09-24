# Serial AMR tree application Mach inlet integration — 2026-09-25

Implementation `4c8844d`, based on `1def375`.

## Change

Fixed and selected serial tree frontends enable Mach input. The driver skips
reservoir construction on active Mach faces, uses the specified inlet temperature,
masks inlet pressure for inactive faces and forwards controls through coupled
common/recursive stepping. Each step constructs current packed tree ghosts and
limits hydro and molecular transport from those states. Subcycling converts
local bounds to the coarse interval using cumulative refinement ratio for hydro
and its square for diffusion. Existing retry and wall bounds remain.

The existing checkpoint context serializes the full configuration, so it already
binds Mach, inlet pressure and temperature without a new checkpoint format.
MPI gas and spray application input remain gated pending their integration.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug. Fixed and selected elementary tree applications
built and executed. `regression_amr_tree_mach_restart_3d` passed: **1/1 tests**,
no failures or skips, 107.68 s (CTest total 107.72 s).

Both mechanisms use Mach 0.2/900 K, extrapolated pressure, reaction and molecular
transport with recursive stepping and three scheduled regrids. Continuous and
checkpoint/resumed final state and budget CSVs match byte for byte; each finishes
three steps at 1e-9 s. Changed Mach, inlet temperature, inlet pressure, physics,
schedule or timestepping rejects restart. Existing truncated checkpoint, invalid
schedule position, inconsistent ledger and output-alias guards also pass.
Provenance, checker Python syntax, project contract and whitespace passed.
No tolerance changed. This is bounded application/restart evidence, not physical
nozzle or ignition validation; the short fixture does not stress the timestep cap.

## Commands and artifacts

Build log: `build/mach-inlet-amr-app-checks/build.log`.
Application logs and checkpoints:
`build/spray-amr-reactive-debug/tests/amr_tree_mach_restart_3d/{fixed,selected}/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_amr_tree_reactive_3d test_selected_elementary_amr_tree_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^regression_amr_tree_mach_restart_3d$' --output-on-failure
python3 -m py_compile tools/check_amr_tree_restart_3d.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Environment overrides are host-local workarounds, not application defaults.
