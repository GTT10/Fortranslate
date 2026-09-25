# Signed PPM transverse face predictors — 2026-09-26

Implementation `4baf6343`, parent `840d399e`, local main, version 0.246.0.
Pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` Godunov.H
`pc_transdo_passive` and `pc_transdd_passive` divide corrected species densities
by corrected total density without clipping negative predictor weights.
PeleF's regular pressure-traced transverse helper rejected such intermediates,
preventing its already signed-aware primitive Riemann path from receiving them.
That predictor-only rejection is removed. Density/finite checks remain.
Accepted conserved cell states still pass the ordinary species/EOS checks before
publication; the timestep transaction is unchanged. No tolerance changed.

## Focused evidence

GNU Fortran 13.3.0, Python 3.12.3, existing Ninja Debug tests-enabled
`build/spray-amr-reactive-debug`. The extended halo unit failed on the old
implementation at `signed transient species after transverse correction`.
After the numerical fix, the initial test incorrectly required a negative
species flux. The pinned Riemann star-state construction clips star species,
so a negative face input does not imply a negative flux. That incorrect test
assertion was replaced with successful Riemann evaluation and mass-flux closure.
An intermediate pressure change did not resolve that assertion; final pressure
is 300000 Pa. All failed logs are retained.

The final unit covers one- and two-direction corrections with weights -0.001
and 1.001, unchanged predictor values, species-sum closure, subsequent Riemann
success/flux closure, and rejection by the ordinary cell-composition EOS entry.
Existing invalid-density rollback and halo/subdomain/physical-boundary checks
remain. Final halo check **1/1 passed**, 6.33 s test / 6.38 s total. The PPM step
check passed separately in the initial combined run, **1/1**, 3.59 s; its solver
source did not change afterward. Thus two distinct relevant tests passed,
without claiming the initial combined run passed. No skips.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_halo_3d test_reactive_ppm_step_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(unit_reactive_ppm_halo_3d|unit_reactive_ppm_step_3d)$' --output-on-failure
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_halo_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, Python syntax, contract and whitespace passed. Artifacts under
`build/ppm-signed-transverse-checks/`: pinned `Godunov.H`, `red.log`, build logs,
initial `ctest.log`, intermediate `ctest-fixed.log`, diagnostic `debug.log`,
final `ctest-final.log` and `unit.log`. Build warnings remain.
This is source alignment and focused red/green evidence, not a compiled full-CTU
comparison, full configuration matrix, application rerun or physical validation.
