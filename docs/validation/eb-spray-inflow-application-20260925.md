# Serial EB spray inflow application — 2026-09-25

## Implementation

Source `ad0e024`, parent `3f78edb`, local main. Both fixed and selected serial
EB spray frontends opt into explicit boundary input. The shared application
accepts only inflow/outflow and requires particle_outflow=true. It constructs
inlet conserved states using the initial mixture, forwards them to the coupled
step, and includes one-cell exterior donors in hydro/diffusion timestep bounds.
Checkpoint context appends face labels and all inlet density/pressure/velocity
values only for explicit boundaries; legacy serialization is unchanged.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug, Python 3.12.3. Existing fixed tests-enabled
`build/spray-amr-reactive-debug` and selected tests-disabled
`build/spray-amr-reactive-selected-debug` configurations; selected bundle
`mechanisms/h2o2_full.json`. No compiler/backend changes or install audit.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_reactive_eb_3d_selected_spray --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_eb_spray.py --exe build/spray-amr-reactive-debug/pelef_eb_spray_3d --source . --work build/eb-spray-inflow-application-checks/fixed --inflow --restart
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/application/check_eb_spray.py --exe build/spray-amr-reactive-selected-debug/pelef_reactive_eb_3d_selected_spray --source . --work build/eb-spray-inflow-application-checks/selected --inflow --restart --selected
python3 -m py_compile tests/application/check_eb_spray.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Each checker passed continuous/save/resume (3 successful executions) and
10 expected rejection executions: changed hydro, inlet velocity, face labels,
unsupported wall, closed particles, checkpoint overwrite, spray context,
geometry context, truncated payload and output alias. Gas/particle/film CSVs
and printed cumulative ledgers match exactly between continuous and resumed
runs. No unexpected runtime failures or skips. Empty film CSV equality does
not qualify active films.

Two additional fixed first-step runs use cfl=0.0001, stop_after_checkpoint=true
and inlet y-speed 200/400 m/s. Accepted dt decreased from
4.012848354806156e-11 to 3.714716476759447e-11 s, below maximum_dt=1e-8.
Inputs and logs are retained as bound_200/bound_400; results in
boundary_timestep.json. These runs check the current-donor bound is active.

Provenance, syntax, contract and diff checks passed. Mechanism SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
The first fixed build invocation used a nonexistent target and failed before
compilation; the corrected target above passed. Compiler warnings remain in
existing code. No test tolerances or golden data were changed.

Artifacts under `build/eb-spray-inflow-application-checks/`:

- build-fixed.log, build-selected.log and initial-target-error.log.
- fixed/eb_spray_glojjap4: inputs, logs, CSVs, checkpoint and timestep comparison.
- selected/eb_spray_2_i8pefx: inputs, logs, CSVs and checkpoint.

## Remaining scope

No periodic/domain-wall particle dispatch, independent inlet mixture,
MPI/tree boundary application dispatch or active film qualification is added.
No full configuration matrix, external solver parity or physical nozzle/spray
validation was run. These focused application checks do not establish complete
PeleC parity.
