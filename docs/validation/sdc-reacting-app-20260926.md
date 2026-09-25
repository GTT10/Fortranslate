# Reacting SDC application routing — 2026-09-26

Implementation `6f445d74`, parent `4fb2df38`, local main, version 0.246.0.

Fixed and selected regular gas applications now allow periodic PPM/PeleC SDC
with chemistry enabled. The shared driver retains reaction-source estimates
between time steps and passes configured chemistry tolerances to the forced
implicit reactor. Nonreacting runs forward absent reaction arguments, retaining
their previous route. Molecular transport and external force can participate
at the same time. Selected-mechanism diagnostic output identifies the actual
forced adaptive implicit SDC integrator. Example:
`cases/reactive_boundaries_3d/sdc_reacting.nml`.

## Verification scope

Completed: both builds passed; **4 matched-composition application cases
passed**. Chemistry changes temperature by **0.004657404120280262 K** in both
variants. Maximum scaled difference across all reacting fixed/selected CSV
entries is **0**. The normalized Euler inventory residual is **1.66121e-17**.
Provenance, Python syntax, project contract and whitespace passed. No final
application failures or skips; existing build warnings remain.

GNU Fortran 13.3.0 Debug builds use fixed tests-enabled
`build/spray-amr-reactive-debug` and selected tests-disabled
`build/spray-amr-reactive-selected-debug`, with the pinned H2/O2 bundle.
The application comparison uses a uniform 4³ mixture at 1100 K and 300000 Pa,
transport enabled, external force (20000,-30000,40000) N/m³, two SDC iterations,
and six time steps covering 1e-8 s. It compares chemistry enabled/disabled
and fixed/selected output CSVs at matching compositions.

The initial cross-variant comparison failed because fixed inputs resolve
`x_h2`, `x_h`, etc., whereas selected inputs use named composition arrays.
Initial selected outputs are retained in `initial-composition/`; the selected
fixture was corrected to match the fixed default mixture. No solver tolerance
was changed. Initial runs are not counted as matched-composition evidence.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_3d_selected --parallel 2
python3 build/sdc-reacting-app-checks/check.py
python3 -m py_compile tools/check_project_contract.py build/sdc-reacting-app-checks/check.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/sdc-reacting-app-checks/` contains checker, inputs, CSVs,
build/run logs, provenance and comparison report. No CTest rerun was required
for this application dispatch change; library SDC tests are recorded separately
in `sdc-reactions-20260926.md`. Uniform transport does not establish nonuniform
reaction/diffusion accuracy. Physical boundaries, AMR/MPI/EB SDC and external
solver/physical validation remain open. These applications have no restart
interface, so persistent reaction sources are in-memory state only.
