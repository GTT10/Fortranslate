# Coupled EB CTU inflow/outflow boundaries — 2026-09-25

Implementation `e2eb060`, based on `f564ea2` on `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

The single-level Strang and full reaction/transport drivers now accept explicit
CTU hydro boundary labels and prescribed inflow U/T. Explicit labels require
compatible extended CTU geometry and the `pelec` solver. The full transport path
also requires a caller-supplied stage donor provider; it invokes this provider
four times, twice per transport half interval. Providers must prepare transport
donors consistent with the hydro labels and inflow data.

Actual time-integrated transport fluxes contribute to `boundary_change`, with
apertures and face areas applied once. Hydro wall and chemistry ledgers retain
their distinct meanings. State, temperature and all ledgers commit together.
A failure in the last transport stage restores the complete incoming state and
leaves caller ledger values untouched. The transport-disabled route forwards
boundary data to the Strang hydro/chemistry operation.

Explicit transport boundary labels currently support prescribed inflow/outflow.
Periodic seams and physical wall traction/heat-flux policies remain pending and
reject before stepping. Hydro-only wall and periodic capabilities are unchanged.
Default calls retain outflow hydro and closed transport faces. Application input,
concrete production donor providers, tree/MPI dispatch and full reference parity
remain separate work; this is a library composition increment.

## Focused checks

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug in
`build/spray-amr-reactive-debug`, existing OpenMPI wrapper configuration.
The final targeted run passed **2/2**, no failures or skips, **6.46 s total**:

- `unit_reactive_eb_transport_3d`: 0.99 s, existing transport/full-driver checks.
- `unit_eb_transport_halo_3d`: 5.42 s, existing stage checks plus full-fluid and
  cut-cell CTU composition, all four provider calls, disabled/enabled chemistry
  with resolved chemical exchange, complete boundary/wall/chemistry inventory,
  transport-disabled Strang forwarding, failure at the fourth provider call,
  missing-provider rejection and unsupported wall-transport rejection.

The new coupled inventory check uses `2e-12` scaled by initial component inventory.
Existing tolerances, chemistry bounds and mechanism data were unchanged.
Earlier focused runs also passed before adding the Strang/guard assertions; the
counts above report only the final source run. Compiler warnings remain, including
the existing executable-stack warning for the test's internal provider procedure.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_transport_halo_3d test_reactive_eb_transport_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_(eb_transport_halo_3d|reactive_eb_transport_3d)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, contract and whitespace checks passed. Mechanism SHA:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Artifacts: `build/eb-coupled-boundary-checks/{build.log,build-final.log,LastTest-initial.log,LastTest.log}`.
No full matrix, application run, installed-binary audit, external-solver comparison
or physical validation was performed for this increment.
