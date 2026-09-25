# Conservative periodic EB transport — 2026-09-25

Implementation `fb5a777`, based on `229e1d0` on `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

Single-level EB transport now supports periodic labels through the built-in
SSPRK2 donor path. Velocity derivatives wrap across periodic directions, including
the tangential derivatives on interior and seam faces. Each periodic seam flux is
computed once from the last/first physical cells. Its opposite-face value carries
the same aperture-weighted transfer. Both representations use the minimum limiter
coefficient from both adjacent cells, including species enthalpy flux adjustment.
Periodic aperture pairs must agree to relative roundoff; incompatible pairs reject.

This increment is the transport operator. Periodic configured geometry,
full-gas driver dispatch and application input/restart checks remain pending;
the coupled driver and application still reject periodic input.

## Focused evidence

GNU Fortran 13.3.0, Python 3.12.3, existing Ninja Debug/OpenMPI configuration in
`build/spray-amr-reactive-debug`.

- Initial focused run passed **2/2**, 16.08 s total: existing
  `unit_reactive_eb_transport_3d` 0.97 s and transport halo suite 15.06 s.
- After strengthening the rejection fixture, final `unit_eb_transport_halo_3d`
  passed **1/1**, 14.80 s (14.86 s total), no runtime failures or skips.
- A full-fluid periodic field varying density, pressure, velocity and composition
  matches the regular Fortran transport flux at every face at `3e-12` scaled
  tolerance. The cut-cell case wraps tangentially to the cut plane and retains
  nonperiodic normal faces.
- A deliberately large RHS interval activates species limiting: paired seam
  fluxes remain identical and every active-cell species stays nonnegative within
  `1e-12`. This is an RHS budget test, not a stable large-timestep flow update.
- Bounded SSPRK2 updates close inventory at `2e-12` scaled tolerance, accounting
  for the remaining nonperiodic faces in the cut-cell case.
- A valid x-cut geometry with incompatible periodic x faces rejects without
  changing state, temperature or the previously published integrated flux. The
  final fixture explicitly checks geometry validity so this exercises seam
  compatibility rather than a general geometry failure.

The initial compile found an unintended extra logical array index introduced
while forwarding periodic flags; corrected before tests. Compiler warnings remain,
including the existing test-provider executable-stack warning. No tolerances or
physical-state checks were loosened.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_transport_halo_3d test_reactive_eb_transport_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_(eb_transport_halo_3d|reactive_eb_transport_3d)$' --output-on-failure
# Final rejection-fixture check:
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_transport_halo_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_transport_halo_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, project contract and whitespace checks passed.
Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Artifacts: `build/eb-periodic-transport-checks/{build.log,build-final.log,build-rejection-check.log,LastTest-initial.log,LastTest-final.log}`.
No application, selected-backend matrix, installed-binary audit, external-PeleC
comparison or physical validation was performed for this increment.
