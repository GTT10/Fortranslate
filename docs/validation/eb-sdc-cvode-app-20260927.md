# Single-grid EB CVODE SDC and restart — 2026-09-27

## Changes and identity

Numerical backend repair: `2eb3f538ba514ee6254b5516ad3c3653637bbab4`.
Application integration: `d0e4034a2dfc62a18810c7661e1013355ce52073`.
Parent before these changes: `88b257de0732bcd687e4324d20b674bc4cd9db35`.
Tests ran on the matching combined working tree before the two separate
commits; they are not independent clean-checkout or installed-binary audits.

Fixed and selected `pelef_reactive_eb_3d` frontends now accept
`sdc_reactor='cvode'` with SUNDIALS, reacting chemistry and SDC. `native`
remains the default. Unknown selectors and unavailable backend selection are
rejected by input validation. The application checks callback/selector
agreement and forwards the forced reactor to the EB driver.

The existing checkpoint SDC block uses `SDC_CVODE_SOURCES` for CVODE and
retains `SDC_SOURCES` for native. Reading requires the marker selected by the
current configuration. Native serialization is otherwise unchanged; previous
readers reject the new CVODE marker. No old-checkpoint hash was replaced.

## Numerical defects exposed by integration

The original cut-cell fixture rejected CVODE while native completed. A
second-iteration predictor had density `0.35248998271315335` and species sum
`0.352489982744202`, a relative discrepancy `8.808381844254385e-11`, exceeding
the existing `5e-11` state closure guard.

The forced backend now evolves N-1 species densities. It reconstructs the
largest initial species from the analytic density
`sum(rhoY_initial) + t*sum(external_species)`, following the existing
constant-volume CVODE dependent-species approach. Internal energy still follows
the prescribed source exactly. A one-species mechanism has an analytic density
path with final thermodynamic validation.

Selected H2/O2 then exposed a returned trace density around `-2.204e-39` and
recoverable RHS failures from negative Newton trial coordinates. As in the
existing constant-volume backend, the RHS is extended to nonnegative trial
coordinates before reconstructing the dependent density. This extension does
not publish trial states. At the endpoint, negative independent densities are
rejected below `-16*epsilon*final_density`; machine-scale noise within that
bound is set to zero, and the dependent density restores total mass. Final
nonnegative and thermodynamic checks still apply. This bound is independent of
user integration tolerances; no existing tolerance or golden hash was loosened.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, SUNDIALS 7.2.0 official static Fortran
interfaces. Existing `build/amr-cvode-app-debug` configuration: Debug,
MPI/spray/SUNDIALS ON, selected `mechanisms/h2o2_full.json`. No compiler/backend
configuration change was made. Build warnings include existing real comparisons,
conversions and internal-procedure executable-stack warnings.

| Test | Result |
| --- | --- |
| `application_eb_sdc_cvode` | Passed, 69.33 s |
| `unit_sdc_cvode_bridge` | Passed, 0.19 s |
| `unit_sundials_forced_reactor` | Passed, 0.13 s |

Final run: **3/3 passed**, zero failures/skips, total **69.70 s**.

The application fixture is the existing 4³ plane case, with covered/cut cells,
state redistribution, molecular transport and external forcing. Fixed uses the
elementary mechanism; selected uses full H2/O2 with named H2/O2/N2 composition.
Each checks active species conversion, exact continuous/restarted CSV bytes,
step-1 checkpoint stop, both native/CVODE restart mismatch directions, and
unknown-selector rejection. Rejected input/restart publishes no final CSV.
Normalized inventory residuals (scale is the larger of initial component
inventory and initial mass) were `3.8565617038930043e-15` fixed and
`8.9212520582436895e-16` selected, below `3e-12`.

The reactor unit adds application-tolerance density closure within four machine
epsilons, initially absent species, and the one-species analytic path. Existing
inert-source, energy, failure/reuse and native-comparison checks pass; the early
native/CVODE relative difference is `5.2984497802368217e-11` (criterion `2e-7`).
The bridge test retains source/history and whole-step failure checks.

Initial app attempts failed on the defects above; those failures are preserved
in `initial-failure.log`, `selected-initial-failure.log`, and
`selected-rhs-failure.log`. Temporary diagnostic source edits were restored
before the final build. No earlier failed attempt is counted as passed.

## Commands and artifacts

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/amr-cvode-app-debug --target pelef_reactive_eb_3d pelef_reactive_eb_3d_selected test_sundials_forced_reactor test_sdc_cvode_bridge --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/amr-cvode-app-debug -R '^(application_eb_sdc_cvode|unit_sundials_forced_reactor|unit_sdc_cvode_bridge)$' --output-on-failure
python3 -m py_compile tests/application/check_eb_sdc_cvode_app.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

All final checks completed successfully. The final backend build is in
`build/eb-cvode-app-checks/trial-build.log`; the last unit assertion addition
was built in `unit-build.log`. Other logs and preserved `LastTest.log` are in
`build/eb-cvode-app-checks/`. Final inputs, outputs and checkpoints:
`build/amr-cvode-app-debug/eb-sdc-cvode/eb_cvode_5rhrtehz/`.
Pinned mechanism source SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

This covers single-grid EB applications. Tree/particle input routing, MPI EB
CVODE execution, arbitrary mechanisms, old installed binaries and physical
validation are not established here. AMR application tests predate the backend
repair and were not rerun in this focused set.
