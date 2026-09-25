# Redistributed EB transport halo exchange — 2026-09-25

Implementation `c695db4`, based on `3c0f7b1` on `main`.
PeleC reference remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

The single-level transport RHS and redistributed Euler update now forward
caller-prepared one-cell U/T halos and six logical face switches. They optionally
return time-integrated, species-limited flux per open face area. Apertures and
face areas enter inventory accounting once. The Euler entry publishes its flux
only after successful redistribution and thermodynamic recovery. Failed halo or
redistribution evaluation retains the input state/temperature and the previously
published flux. A successful zero-interval update publishes zero exchange.

This is a stage interface, not physical boundary selection. Callers must rebuild
donors per stage. SSPRK2 donor preparation, periodic seam synchronization, thermal
wall conditions and coupled application dispatch remain pending. Existing calls
without halo arguments retain the closed transport boundary convention.

## Focused evidence

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug in the existing
`build/spray-amr-reactive-debug`; configured OpenMPI 4.1.6 wrappers. No compiler,
MPI or optional-backend configuration changed.

Two tests passed, no failures or skips, **1.48 s total**:

- `unit_reactive_eb_transport_3d`: 0.98 s, existing transport/coupled regression.
- `unit_eb_transport_halo_3d`: 0.46 s, full-fluid and planar cut-cell cases,
  active faces on all axes, independent directional flux comparison, aperture
  inventory after redistribution, covered storage, missing-halo and invalid
  redistribution rollback, and zero-interval exchange. Existing species limiter
  checks also passed. Inventory tolerance remains `2e-12` scaled by initial
  component inventory; no tolerance was loosened.

Commands completed:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_transport_halo_3d test_reactive_eb_transport_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_(eb_transport_halo_3d|reactive_eb_transport_3d)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, Python syntax, project contract and whitespace checks passed.
The implementation differs from the tested source only by explanatory comments.
Mechanism SHA remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Artifacts: `build/eb-transport-redist-halo-checks/{build.log,LastTest.log}`.
This is focused library evidence; no full matrix, application run, external-solver
comparison or physical validation was performed for this increment.
