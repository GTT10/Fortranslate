# Serial tree isothermal transport — 2026-09-26

Implementation 72304fd, parent a6f1a71, local main, version 0.246.0.
Serial common/recursive tree transport accepts wall_temperature(6), validates
finite positive supplied temperatures and forwards them to the patch operator.
Only domain faces receive the physical isothermal label; interior patch faces
retain donor transport. Serial halo reflection treats isothermal walls like
no-slip walls for momentum, retaining donor temperature. Thermal exchange is
computed by the existing wall heat-flux law. Flux synchronization now accepts
the isothermal label. Inputs without wall temperatures preserve their paths.

## Completed checks

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug, existing fixed tests-enabled
build/spray-amr-reactive-debug; no backend changes.
Halo/common-transport/recursive-transport tests passed **3/3**, 0.26, 2.42 and
21.37 seconds (CTest total 24.09 seconds), no failures or skips.
Three plane orientations and four/five-patch trees exercise all-isothermal
exterior walls, 100 K hotter than the hottest initial patch temperature.
Boundary mass/species flux is exactly zero and boundary energy exchange is
positive. Every component inventory residual meets the existing 3e-14 bound;
maximum printed isothermal residual is 6.1811175941805725e-15. Existing
adiabatic, periodic, donor and rollback cases remain in the same tests.
The former unsupported-isothermal halo rejection now uses an unknown wall label,
since isothermal reflection is supported.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_transport_3d test_eb_tree_subcycle_transport_3d test_eb_tree_halo_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_(transport|subcycle_transport|halo)_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, project contract and whitespace passed. Existing compiler
warnings remain. No tolerance or golden-data changes.
Artifacts: build/eb-tree-isothermal-transport-checks/build.log and LastTest.log.
Mechanism SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
Pinned PeleC: bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8.
Pinned PelePhysics: 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.

## Limits

MPI label/temperature propagation, thermal timestep bounds, full reacting split
and application input/restart remain to be connected. These checks use a fixed
small timestep. No full matrix, installed audit, external solver parity or
physical heat-transfer validation ran. Complete PeleC parity remains unproven.
