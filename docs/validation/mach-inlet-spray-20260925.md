# Mach inlet regular spray coupling and restart — 2026-09-25

Implementation `6ab6937`, based on `83a55f3`.
`advance_spray_reactive_3d` appends optional Mach/pressure controls and forwards
them through gas splitting and both static/Dynamic LES half steps. Gas,
particles and optional film/component inventories retain their existing private
candidate transaction. The regular fixed/selected spray namelist accepts
`inflow_mach(6)`: negative retains a reservoir, 0 through 100 requires an inlet
face with zero prescribed velocity. Temperature must be positive; zero pressure
extrapolates the interior, positive pressure up to 1e15 Pa prescribes it.

The application skips fixed-state construction for Mach faces, masks reservoir
pressure from the low-level Mach array, and reconstructs padded current inlet
states before each timestep to bound hydro and molecular transport. Existing
coupled retries retain LES stability rejection. A new context marker and the
Mach array bind restart; existing context already binds inlet temperature and
pressure. Disabled Mach inputs retain the prior context encoding. MPI/AMR
Mach dispatch is still open.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug. Fixed and selected H2O2 spray applications built;
selected execution was not run. `spray_dynamic_mach_restart`: **1/1 passed**,
43.27 seconds, total 43.32; no failures or skips.

The 4-cubed fixture uses two evaporating water parcels, molecular transport,
Mach 0.2 at x-low with interior pressure, and Dynamic LES with Gaussian test
filter ratio 4/grid ratio 2. Chemistry is disabled in this application check.
Continuous and stop-after-two-steps/restarted gas and parcel CSVs are byte-exact;
final history ledgers match. Final time is 6e-7 s, with two live parcels.
Mass error is 2.2566521086212286e-23 kg; energy error is
2.8189256484623115e-18 J. Boundary mass and energy are positive and conservation
meets the existing 1e-9 scaled checker bounds. Gas H2 remains exactly zero despite
a supplied H2 reservoir composition, confirming the Mach inlet uses interior
composition. Changed Mach, inlet pressure or temperature rejects checkpoint
continuation before history output; Mach 101 rejects input. Existing Dynamic
model/filter context checks also run in this checker.

Mechanism provenance, Python syntax of both touched checker and contract checker,
project contract and whitespace passed. No tolerances changed. This is regular
spray integration/restart evidence, not selected execution, MPI/AMR integration,
reacting-nozzle parity, physical spray validation or turbulence qualification.

## Commands and artifacts

Logs: `build/mach-inlet-spray-checks/`.
Application artifacts:
`build/spray-amr-reactive-debug/tests/spray_dynamic_mach_restart/spray-boundaries-nm2rjnnu/`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_3d_selected --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_dynamic_mach_restart$' --output-on-failure
python3 -m py_compile tools/check_spray_boundaries.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The MPI-related environment settings are host-local workarounds, not application
defaults.
