# Serial EB tree adiabatic wall halos — 2026-09-26

Implementation eba5239, parent d5b6898, local main, version 0.246.0.
The serial tree halo API now accepts slip_wall and no_slip_wall labels.
One-layer face ghosts reflect adjacent donor momentum. Full/wide halos map
indices through domain-wall reflections before selecting peer/ancestor donors;
slip walls reverse normal momentum, while no-slip walls reverse all momentum.
Repeated reflections accumulate signs, including mixed-wall corners and halos
wider than the domain. Density, energy, species and donor temperature are
preserved. Existing snapshot interpolation and atomic output commit remain.
Periodic and inlet behavior remains covered by the existing test.

## Completed evidence

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug, existing fixed tests-enabled
build/spray-amr-reactive-debug. No backend configuration changes.
The focused unit test passed **1/1**, 0.25 s (CTest total 0.30 s), no failures
or skips. It checks nonzero three-component momentum against analytical sign
and mirrored-index expectations for slip/no-slip faces and six-layer corners.
Unsupported isothermal-wall input leaves state, temperature and volume outputs
unchanged. This does not qualify a thermal wall operator.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_halo_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_halo_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, syntax, project contract and final whitespace checks passed.
A trailing space in the test was removed after the first diff check.
Compiler warnings remain; no numerical tolerance or golden-data changes.
Artifacts: build/eb-tree-wall-halo-checks/build.log and LastTest.log.
Mechanism SHA: 0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4.
Pinned PeleC: bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8.
Pinned PelePhysics: 6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b.

## Remaining work

This qualifies serial halo donor construction only. Reflected extended EB
geometry, wall-aware coupled stepping, MPI wall halos, application controls
and isothermal tree walls require separate connection and checks. No full
matrix, installed audit, external solver parity or physical validation ran.
Complete PeleC parity remains unproven.
