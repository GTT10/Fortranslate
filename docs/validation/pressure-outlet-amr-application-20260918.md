# Serial two-level AMR pressure outlet application — 2026-09-18

Regrid implementation `23e5f32`, application/restart `6ecfe2b`, base `887ee83`.

Fixed and selected serial two-level AMR front ends enable pressure input and
forward targets into the complete reactive split. Regridding now accepts the
same targets and fills coarse interpolation ghosts with the pressure kernel.
The application supplies physical faces for this pressure-aware interpolation;
the existing tagging policy remains unchanged.

A `PRESSURE_OUTLET_V1` suffix and all six targets are appended to checkpoint
boundary identity when any target is positive. Existing zero-target identities
remain unchanged. Restoring with a changed or removed target fails before
output publication. Patch trees and MPI applications remain unconnected and
continue to reject nonzero pressure input.

## Focused checks

GNU Fortran 13.3.0 / Open MPI 4.1.6, Debug/Ninja; existing fixed tests-enabled
and selected H2O2 tests-disabled configurations. Selected executable built;
no selected or MPI runtime qualification in this increment.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_amr_reactive_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_amr_reactive_3d_selected --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^amr_pressure_outlet_restart_3d$' --output-on-failure
python3 -m py_compile tools/check_amr_physical_boundary_restart.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The fixture uses molecular transport and inert hydro, a moving isothermal
x-low wall, 102325 Pa x-high pressure outlet, ratio-two boundary-touching patch,
and temperature-triggered regridding. The full run takes more than one coarse
step; a second run stops at the first checkpoint and resumes. Coarse and fine
CSVs must match byte-for-byte, and conserved boundary ledger text must match.
Changing the target to 103325 Pa, removing it, or changing wall temperature
must reject without publishing CSV output. Six executable launches per check.

Artifacts: `build/pressure-outlet-amr-app-checks/`, including provenance,
fixed/selected initial/final build logs, and `ctest*.log`. This is focused
application/restart evidence, not physical or MPI validation.

Final regridding CTest **1/1 passed**, 0 skips, 80.86 s. Both final builds,
provenance, Python syntax, project contract and diff checks passed. Earlier
static-topology fixture passed 1/1 in 71.18 s before the regrid extension.
The first regridding fixture was rejected at input validation because it lacked
`regrid_physical_boundaries`; enabling the required existing policy fixed the
fixture. No solver tolerance or validation rule was loosened.

Final fixture artifacts:
`build/spray-amr-reactive-debug/tests/amr_pressure_outlet_restart_3d/run-glvsw1uc/`.
