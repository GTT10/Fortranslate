# Dynamic LES harmonic interpolation — 2026-09-18

Implementation `e01caa12f648597ccb0e376b4c98ca75f5b0a4a4`, based on `80231c4`.
Serial/MPI dynamic RHS and SSPRK2 APIs accept trailing optional logical
`harmonic_mean`; omission retains arithmetic interpolation. The option applies
to dynamic coefficients, deviatoric stress quantities, isotropic strain quantity
and thermal flux quantity. Velocity, temperature and composition retain arithmetic
face averages. MPI checks option agreement before coefficient/halo communication,
including when dt is zero. Coupled-spray/application selection remains pending.

Reference: pinned PeleC `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`,
`PeleCUtilities.H::pc_move_transcoefs_to_ec`, used by dynamic LES in `LES.cpp`.
Source SHA256 is in `references/les_face_interpolation.json`; downloaded source
is `build/les-reference/PeleCUtilities.H`. Same-sign values have signed harmonic
mean; opposite signs or either zero produce zero. The elemental helper uses
min/max magnitude scaling, avoiding the reference formula's product overflow
and underflow. Very small same-sign inputs therefore retain their mathematical
mean instead of becoming zero through product underflow. This is a deliberate
numerical difference, not bitwise C++ parity.

## Independently completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Ninja Debug; unchanged configurations.
Fixed tests-enabled `build/spray-amr-reactive-debug` built `test_spray_core` and
`test_les_mpi_dynamic`. Selected H2O2 tests-disabled
`build/spray-amr-reactive-selected-debug` built `pelef_selected_spray_runtime`.

Focused CTest **4/4 passed**, 1.45 seconds, zero failures/skips. Scalar checks
cover positive/negative harmonic means, opposite signs, zero, arithmetic mean,
and magnitudes `1e300`/`1e-300`. Two/four-rank tests exercise arithmetic and
harmonic variants for box, Gaussian, box5opt and gaussian3opt: exact serial RHS,
diffusion estimate, state and temperature parity; conservation; dt rejection;
and rank-dependent interpolation refusal with exact RHS preservation. Harmonic
flux differs from arithmetic on the nonuniform fixture. The existing serial
SSPRK2 test also passes. Provenance, Python syntax for the contract checker,
project contract and whitespace checks passed. No failed builds/tests occurred.
Selected runtime was compiled, not exercised. No full matrix, external solver
parity or physical validation ran.

## Commands and artifacts

Logs: `build/les-harmonic-checks/{provenance,build,build-selected,ctest}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core test_les_mpi_dynamic --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_les_(face_mean|dynamic_step|mpi_dynamic_[24])$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
