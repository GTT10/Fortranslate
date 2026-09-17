# Serial EB impact application and restart — 2026-09-18

Implementation `266f495`, base `10280b6`.

Fixed/selected serial EB spray applications now accept `wall_impacts` and
`spray_eb_wall_options` through `walls` in the namelist. The allocated reactive
path retains growing/consumed populations and evolving stationary films. A
separate film CSV records geometry, temperature and total component masses.
The existing checkpoint format retains film/component arrays and uses its
retired-parent mode for consumed splash parents. Restart validates ordered IDs,
watermark/injection/splash consistency, fluid particle locations and EB film
geometry. Context V2 binds all settings and rejects V1 snapshots. Film-output
paths participate in input and checkpoint alias rejection.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja, existing fixed tests-enabled and
selected H2O2 tests-disabled build directories.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_(app_deposit|app_splash|restart)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_eb_spray.py --exe build/spray-amr-reactive-selected-debug/pelef_reactive_eb_3d_selected_spray --source . --work build/spray-eb-impact-app-checks/selected --selected --scheduled --restart --impacts splash
python3 -m py_compile tests/application/check_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final fixed CTest **3/3 passed**, total **4.53 s**. Selected splash script passed
in `selected/eb_spray_14pka4cm/`. Provenance, syntax, contract and whitespace
passed. Gas/particle/film CSV bytes and printed cumulative ledgers match exactly
between uninterrupted execution and continuation. Impact snapshots are saved
after step two, after impact processing; the legacy specular fixture keeps its
step-one split. Existing overwrite, changed input/geometry, truncation and alias
rejections also run. No skips.

Five scheduled events include tied and zero-dose events and final-time injection.
The deposition fixture uses 20 m/s and template mass 1e-10 kg. The final splash
fixture uses 200 m/s and 1e-16 kg to keep repeated colocated impacts within the
thin-film correlation range; it retains 24 children and two final-time injected
parents with final IDs 31--32. Deposition retains two final-time injected parents
and positive film inventory. Source substeps are bounded to eight in the impact
fixtures. These are small synthetic interface/restart checks, not physical spray
or general wet-impact validation.

## Failures and remaining limitation

The initial splash fixture used 1e-10 kg and timed out at 60.11 s; the other two
tests passed. A debugger reproduction confirmed that the existing wet splash
mass partition at alpha=pi/2, film ratio=0.11228532717131677 and splashed=1 gives
fractions [0.252,3.3463,0.252,-2.8503]. The kernel correctly rejects this negative
mass fraction. The initial application repeatedly reduced dt toward the same
inadmissible impact. The correlation and conservation checks were not changed;
this thick-film condition remains unsupported. The debugger reproduction was
stopped at the invalid partition and is not a completed simulation.

The first thin-film run exposed checkpoint rejection because consumed parents
were absent. Enabling the checkpoint API's existing retired-parent mode fixed
the writer/reader contract. This was independently rerun successfully.

Artifacts: `build/spray-eb-impact-app-checks/` contains final and initial CTest
logs, checkpoint-failure log, build logs, selected run, provenance/contract logs,
and `wet-diagnostic.log` with the invalid fractions. No full matrix, installed
audit, external parity, convergence or physical validation was run. MPI impact
and film ownership remain open.
