# Serial EB spray open-boundary coupling — 2026-09-25

## Scope and source

Implementation: `6e70b92`, based on `ffe3f4f`, on local `main`.
`advance_eb_spray_reactive_allocated` forwards explicit inflow/outflow faces
and inlet conserved states/temperature guesses into the existing CTU gas split.
Explicit faces require extended geometry and enable particle outflow on all
open faces. Gas inflow does not create particles. An explicitly false
`particle_outflow` conflicts with these faces and is rejected. Periodic and
physical wall labels are rejected pending consistent particle dispatch.
The fixed-array wrapper and application inputs are unchanged.

Both particle halves contribute escaped inventory to `boundary_change`.
The existing drift path retains escaped parcels with zero mass and clears
component masses. State, particles and ledgers publish only after success.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug, Python 3.12.3; existing tests-enabled
`build/spray-amr-reactive-debug` configuration. No compiler or backend change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_reactive(_ctu)?$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest: **2/2 passed, 0 failed, 0 skipped**, 2.41 s total
(legacy 0.30 s, CTU 2.05 s). The CTU fixture checks coupled chemistry,
transport, evaporation and EB reflection with an open gas inlet; total
inventory closure uses the existing 3e-12 relative tolerance. It also checks
unsupported-label rollback, invalid inlet density after the first particle
half, and inventory loss when a parcel escapes through the inlet.
Provenance, Python syntax, project contract and whitespace checks passed.
The pinned mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

Two initial test attempts failed because of fixture assumptions: a negative
inlet temperature guess is recoverable from the conserved state, and escaped
parcels retain zero-mass slots instead of being removed. The fixture now uses
negative density and checks zero represented mass. No numerical tolerances
were changed. Existing real-comparison compiler warnings remain.

Artifacts: `build/eb-spray-open-checks/` contains `initial-build.log`,
`build.log`, both failed-attempt logs and final `LastTest.log`.

## Limits

This is a fixed-mechanism serial library check. Application input/restart,
selected-mechanism execution, explicit-boundary impact/film coupling,
periodic/wall particle dispatch and MPI/tree dispatch are not qualified here.
No full matrix, installed-binary audit, external solver comparison or physical
spray validation was run. Full PeleC port completion is not established.
