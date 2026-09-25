# Isothermal EB patch transport — 2026-09-26

Implementation b02b1a1, parent 4959897, local main, version 0.246.0.
The transactional patch Euler transport API now accepts optional wall_temperature
and forwards it with physical-face labels to existing isothermal wall fluxes.
This is the lower-level connection needed by tree transport. Existing callers
without wall temperature retain their behavior.

## Completed check

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug and existing OpenMPI wrappers.
Fixed tests-enabled build/spray-amr-reactive-debug; no backend changes.
The expanded transport-halo test passed **1/1**, 16.38 s (CTest total 16.43 s),
no failures or skips. Regular and cut-cell geometries exercise each open
orientation/side with a moving gas state and wall 100 K hotter than gas.
Patch integrated flux agrees with dt times the public flux operator within
3e-12 relative scale. Gas energy increases. Negative wall temperature rejects
without changing state, temperature, limiter or any integrated face flux.
The same test retains its existing regular-transport wall-flux reference and
impermeability checks. No tolerance or golden-data changes were made.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_transport_halo_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_transport_halo_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, project contract and whitespace passed. Compiler warnings
and the test object's executable-stack linker warning remain.
Artifacts: build/eb-isothermal-patch-checks/build.log and LastTest.log.
Mechanism SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
Pinned PeleC: bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8.
Pinned PelePhysics: 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.

## Limits

Tree/MPI wall-temperature forwarding, thermal timestep restrictions and tree
application input/restart remain to be connected. No full matrix, installed
audit, external solver parity or physical validation ran. Complete PeleC
parity remains unproven.
