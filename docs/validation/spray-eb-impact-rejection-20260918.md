# EB impact correlation rejection — 2026-09-18

Implementation `f8977a8`, base `30027f8`.

The negative thick-film mass fraction is present in the pinned upstream
[AhamedSplash.H](https://github.com/Pele-Suite/PelePhysics/blob/6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b/Source/Spray/BreakupSplash/AhamedSplash.H#L88-L116).
The local reference copy's SHA-256 was checked again:
`37ddd95a06bc96d6c2b1e740cf583a872335efc787b56a4cce80211f4381cb92`, matching
`references/spray_wall_impact.json`. For normal incidence the thick-film branch
returns [0.252,3.3463,0.252,-2.8503] times the splashed fraction. This is not a
Fortran coefficient transcription error. Correlations and admissibility checks
are unchanged; general thick-film physical support remains unresolved.

The EB contact adapter now exposes `model_rejected` for invalid empirical
classification, mass partition or outgoing normal speed. Impact population,
source and reactive composition propagate the flag while preserving state.
The application reports the model-range error immediately instead of reducing
dt repeatedly into the same inadmissible collision. Other step failures retain
their original retry behavior. No alternate impact physics is substituted.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja, existing fixed tests-enabled and
selected H2O2 tests-disabled configurations.

```sh
sha256sum build/wall-impact-reference/AhamedSplash.H
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d test_spray_eb_impact_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_(impact_reactive|app_splash)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_eb_spray.py --exe build/spray-amr-reactive-selected-debug/pelef_reactive_eb_3d_selected_spray --source . --work build/spray-eb-impact-rejection-checks/selected --selected --scheduled --restart --impacts splash
python3 -m py_compile tests/application/check_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Fixed CTest **2/2 passed**, total **3.38 s**. Selected application script passed
in `selected/eb_spray_qdrw46hj/`. Builds, provenance, Python syntax, contract and
whitespace passed. No failures or skips in this increment.

The reactive library test supplies a preexisting film, confirms the rejection
flag and zero accepted steps, and checks gas, component arrays, film mass,
watermark and diagnostics remain unchanged. The application script retains the
thin-film success/restart checks, then restores the formerly timing-out 1e-10 kg
injection mass. Both fixed and selected applications now exit with the specific
model-range message within the negative test's ten-second subprocess limit.
The complete fixed application test, including all its other runs, took 1.83 s.

Artifacts: `build/spray-eb-impact-rejection-checks/` contains build, CTest,
selected, provenance and contract logs. Each application directory includes
`wet_rejection.log`. No broad matrix, external solver parity, installed audit or
physical validation was run. This fixes failure handling, not the correlation's
physical validity domain.
