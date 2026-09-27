# Prescribed periodic transport coefficients, 2026-09-27

Implementation `23174a8918fb7662e49cd7b25f4f82858e19bef8`, based on
`3e0f4957`, local `main`.

## Capability

`reactive_transport_face_flux_3d`, `reactive_transport_fluxes_3d` and
`reactive_transport_timestep_3d` now accept optional
`constant_coefficients=[dynamic_viscosity,thermal_conductivity]` in SI units
(Pa s and W/(m K)). The supplied coefficients replace the mixture values
before stress, viscous work, heat flux and diffusion timestep calculations.
They must be finite and nonnegative. Omission preserves the existing molecular
transport path. Existing enable flags still control the corresponding terms.

This API supplies part of the constant-property path needed by the HIT
reference. Species diffusion coefficients remain the existing mixture model;
disable species diffusion to obtain the reference's zero species diffusion.
Bulk viscosity remains zero. Existing thermodynamic and molecular database
validation is retained, including calculation of face composition/enthalpies.
This is not a new mixture database format and does not alter generated data
or checkpoint identity. Physical-boundary/AMR callbacks, SDC/application input
routing, Reynolds-based coefficient construction and complete HIT parity are
not yet supplied by this increment.

## Checks

Used existing `build/amr-cvode-app-debug`: GNU Fortran 13.3.0 Debug,
MPI/spray/SUNDIALS enabled, Open MPI 4.1.6, SUNDIALS 7.2.0. No compiler or
backend configuration changed. The focused test links the fixed mechanism.
Existing real-comparison and executable-stack warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/amr-cvode-app-debug \
  --target test_constant_transport_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/amr-cvode-app-debug \
  -R '^unit_constant_transport_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Provenance, build, syntax, contract and whitespace checks pass.
CTest: **1/1 passed, 0 failures, 0 skips**, 0.05 s (total 0.10 s).

- All three face normals match independent Newtonian deviatoric stress,
  viscous work and Fourier heat-flux formulas; maximum absolute flux error is
  **0**, against `1e-12`.
- A nonuniform periodic 4³ state verifies exact doubling of all three flux
  arrays when both coefficients double. Maximum diffusivity doubles and the
  stable timestep halves exactly.
- Zero coefficients yield exactly zero stress/heat flux with species diffusion
  disabled. Negative flux coefficients and nonfinite timestep coefficients are
  rejected with zero outputs.

The focused test verifies flux and timestep kernels, not an evolved application
conservation run or external solver parity. No tolerances or golden data changed.
Artifacts are `build/constant-transport-checks/{build.log,ctest.log,LastTest.log}`.
The pinned mechanism YAML retains SHA-256
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
