# Serial EB spray domain-wall coupling — 2026-09-26

## Implementation

Source f88e021, parent f3e9ec6, local main. The allocated reacting-spray API
accepts slip_wall, no_slip_wall and isothermal_wall gas faces, forwards wall
temperatures to the gas split, and applies specular particle reflection at
those same outer faces. Particle mass and energy are unchanged by reflection;
the normal momentum change joins the exterior inventory. Embedded-surface
exchange retains its existing wall ledger. Gas no-slip does not impose
particle tangential sticking. Impact/film options with periodic or domain-wall
faces reject because their event handler supports only outer outflow.

The shared fixed/selected application now accepts these labels, passes thermal
controls, adds the existing EB wall heat-transfer timestep bound, and appends
wall-v1 plus six wall temperatures to isothermal checkpoint context. Other
boundary context remains unchanged.

## Completed verification

GNU Fortran 13.3.0, Ninja Debug, Python 3.12.3. Existing fixed tests-enabled
build/spray-amr-reactive-debug and selected tests-disabled
build/spray-amr-reactive-selected-debug (H2O2 full bundle). No backend change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_reactive(_ctu)?$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_eb_spray.py --exe build/spray-amr-reactive-debug/pelef_eb_spray_3d --source . --work build/eb-spray-wall-coupling-checks/fixed --wall --restart
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_eb_spray.py --exe build/spray-amr-reactive-selected-debug/pelef_reactive_eb_3d_selected_spray --source . --work build/eb-spray-wall-coupling-checks/selected --wall --restart --selected
python3 -m py_compile tests/application/check_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Library tests: **2/2 passed**, 4.42 s total (legacy 0.29 s, CTU 4.09 s).
The CTU case reflects a parcel at each of the three gas wall types with
chemistry, transport and evaporation active. Gas-plus-particle inventory
closes at the existing 3e-12 relative tolerance. Invalid wall temperature
rejects without publishing state, component mass or boundary ledger. Previous
open/periodic/force cases remain in the same test.

Each fixed/selected application check passed **3 successful runs and 9
expected rejections**. At 1200 K y/z walls, injected particles reverse their
y velocity and remain within the domain. Continuous and resumed gas, particle,
empty-film CSVs and cumulative ledgers match exactly. Rejections cover changed
wall temperature, invalid temperature, changed wall type, changed hydro,
checkpoint overwrite, spray context, geometry context, truncated payload and
output alias. No unexpected failures or skips occurred.

Provenance, syntax, project contract and whitespace checks passed. Mechanism
SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
No tolerances or golden data changed. Existing compiler warnings remain.

Artifacts in build/eb-spray-wall-coupling-checks/: build.log, build-fixed.log,
build-selected.log, LastTest.log, fixed/eb_spray_37kkmalw and
selected/eb_spray_8el0grck. The library log predates application-only edits;
its tested numerical wrapper is unchanged in the final commit.

## Limits

No external-wall sticking/splash/film model or MPI/tree dispatch is added.
No full matrix, installed audit, external solver parity or physical spray
validation was run. The application run exercises thermal-bound invocation;
no separate active-wall timestep scaling comparison was rerun here.
Full PeleC parity is not established.
