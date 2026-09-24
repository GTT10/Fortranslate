# Mach inlet regular application input — 2026-09-25

Implementation `41f46a9`, based on `73d716e`.
The regular fixed and selected 3D reacting-flow front ends accept
`inflow_mach(6)` in `reactive_3d`: negative disables the mode, 0 through 100
selects an inward normal Mach inlet on an `inflow` face. Existing
`inflow_temperature` prescribes temperature; `inflow_pressure=0` extrapolates
interior pressure and a positive value prescribes pressure (default 101325 Pa,
maximum 1e15 Pa). Conflicting nonzero `inflow_velocity` is rejected.
Composition is taken from the current stage's interior donor. Shared reservoir
composition continues to apply to ordinary reservoir faces.

The application bypasses fixed-reservoir construction for Mach faces, forwards
controls through the full reactive composition, and reconstructs current padded
inlet states before every step for hydro and molecular diffusion timestep limits.
Only active Mach inputs allocate this additional timestep workspace. Reservoir
faces can coexist; their pressure values are masked out of the Mach control
array. The reader's new explicit capability flag defaults false: AMR/MPI and
other unconnected callers reject active Mach settings. Its existing pressure
outlet capability gate remains independent. No restart format was changed.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug. Fixed application/unit and selected H2O2
application builds passed. Selected execution was not run.
Three focused CTests passed, with no failures or skips:

- `unit_reactive_boundaries_3d`: 7.96 s, including reader capability gates,
  extrapolated/prescribed inlet pressure, velocity conflict, face and Mach-range
  validation, plus existing six-face hydro/ghost checks.
- `reactive_boundary_3d_inflow`: 0.85 s, existing reservoir regression.
- `reactive_boundary_3d_mach_inlet`: 97.21 s, full H2O2 chemistry and all molecular
  processes, PLM/HLLC, x-low Mach 0.2 with extrapolated pressure and z-high Mach
  0.1 at 120000 Pa. Completed two steps to 2e-7 s. Boundary-accounted normalized
  Euler conservation error 1.5720931501039814e-18, species closure error zero,
  minimum transport theta 1. Mass increase 1.2886910957610153e-10 kg and energy
  increase 1.6564629694351858e-4 J agree with boundary exchange within the existing
  application gate.

Total CTest time 106.07 s. Mechanism provenance, Python syntax of the contract
checker, project contract and whitespace passed. No tolerances changed.
This is a small regular application integration check, not external nozzle
parity, physical validation, scaling qualification or LES/AMR/MPI integration.

## Commands and artifacts

Logs: `build/mach-inlet-input-checks/` (fixed/selected build, CTest and copied
`LastTest.log`). CSV:
`build/spray-amr-reactive-debug/tests/boundary_3d_mach_inlet/mach_inlet.csv`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_boundaries_3d pelef_reactive_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(unit_reactive_boundaries_3d|reactive_boundary_3d_(mach_inlet|inflow))$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The MPI-related environment settings are host-local workarounds, not application
defaults.
