# Built-in EB inflow/outflow transport donors — 2026-09-25

Implementation `84c0857`, based on `f514f3a` on `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

The single-level transport SSPRK2 API accepts explicit inflow/outflow labels and
prescribed inlet U/T. It calls the existing regular boundary filler with one-cell
halos at each stage. Direct boundary data and a custom provider are mutually
exclusive in this API. The full gas driver selects the built-in route when labels
are present without a custom provider, using the same inlet data as CTU hydro.
Custom providers remain available. Boundary exchange accounting applies to both
routes. No new host-associated procedure is passed by the production driver.

Defaults remain outflow hydro with closed transport faces. Periodic and physical
wall transport labels still reject pending their flux policies. Missing prescribed
inflow data rejects transactionally. Input-file selection, inlet timestep limits,
restart identity and application/tree/MPI dispatch remain pending.

## Focused evidence

Existing Ninja Debug directory `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, Python 3.12.3, existing OpenMPI wrapper configuration.
Final targeted run: **2/2 passed**, no failures or skips, **8.99 s total**:

- `unit_reactive_eb_transport_3d`: 0.98 s, existing default transport/full-driver checks.
- `unit_eb_transport_halo_3d`: 7.96 s, full-fluid and cut-cell coupled CTU with
  disabled/enabled chemistry. Built-in donors produce exactly the same U, T,
  boundary, wall and chemistry ledgers as the independently supplied stage
  provider. Existing stage inventory and failure rollback checks also pass;
  the former missing-provider rejection now checks missing prescribed inlet data.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_transport_halo_3d test_reactive_eb_transport_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_(eb_transport_halo_3d|reactive_eb_transport_3d)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, Python syntax, contract and whitespace checks passed. The committed
implementation differs from the tested source only by one indentation correction.
Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Artifacts: `build/eb-built-in-donor-checks/{build.log,LastTest.log}`.
Compiler warnings remain, including the existing executable-stack warning for the
test's internal custom provider. No full matrix, application run, installed-binary
audit, external-solver comparison or physical validation was performed here.
