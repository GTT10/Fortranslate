# Serial EB tree adiabatic wall hydro — 2026-09-26

Implementation 88274e3, parent af56f78, local main, version 0.246.0.
Serial tree flux synchronization now accepts slip_wall and no_slip_wall labels.
This connects the reflected halo donors and existing physical CTU wall flux
handling to common and recursive hydro updates without rejecting final reflux.
Periodic flux synchronization is unchanged; wall labels do not create seam peers.

## Completed focused checks

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug, existing fixed tests-enabled
build/spray-amr-reactive-debug. No backend changes.
Hydro and recursive-hydro tests passed **2/2**, 36.57 and 4.35 seconds
(CTest total 40.97 seconds), no failures or skips.

Three plane orientations and four/five-patch trees exercise common PCM and CTU
with mixed tangential slip/no-slip walls, and recursive PCM with tangential
no-slip walls. All component inventory residuals satisfy the existing 3e-14
bound after subtracting boundary and embedded-wall ledgers. Largest printed
recursive-wall residual: 6.3574333142149996e-15. Existing periodic, outflow,
Mach inlet and transactional failure checks remain in those tests.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_hydro_3d test_eb_tree_subcycle_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_(hydro|subcycle_hydro)_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, project contract and whitespace passed. Existing compiler
warnings remain. No tolerances or golden data changed.
Artifacts: build/eb-tree-wall-hydro-checks/build.log and LastTest.log.
Mechanism SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
Pinned PeleC: bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8.
Pinned PelePhysics: 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.

## Limits

Walls are tangential to the planar EB, so existing extended plane geometry
is compatible. General reflected geometry, recursive CTU wall qualification,
MPI wall hydro, viscous/thermal wall transport, full reacting split, application
inputs and isothermal tree walls remain separate work. No full matrix,
installed audit, external solver parity or physical wall-flow validation ran.
Complete PeleC parity remains unproven.
