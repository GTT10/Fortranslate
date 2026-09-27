# Molecular-data-independent prescribed transport, 2026-09-27

Implementation `165bde75b9c02b038898c3a4a67a27c5d327a273`, based on
`ee15d193`, local `main`.

## Change

The explicit viscosity/conductivity path previously calculated molecular
coefficients before replacing their values. This unnecessarily required a
Lennard-Jones database even when species diffusion was disabled. Such a
requirement obstructs connecting a calorically perfect pseudo-gas to HIT.

The periodic 3D face and timestep kernels now skip molecular-coefficient
calculation when explicit coefficients are present and species diffusion is
disabled. Face thermodynamic properties and compositions are still validated;
unused species flux/enthalpy work arrays are zeroed. The transport array must
still have the correct species dimension, but its molecular values are not
used in this path. Enabling species diffusion or omitting the coefficients
retains the existing molecular path and validation. No generated mechanisms,
checkpoint records, bounds or tolerances changed. Gamma-law application
integration remains to be implemented; this change removes its unnecessary
molecular-data prerequisite.

## Checks

Existing `build/amr-cvode-app-debug`, GNU Fortran 13.3.0 Debug, Open MPI 4.1.6,
SUNDIALS 7.2.0, MPI/spray/SUNDIALS enabled. No compiler or optional-backend
configuration changed. The test links the fixed-mechanism core.

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

Final CTest **1/1 passed**, 0 failures, 0 skips, 0.02 s (total 0.07 s).
The existing analytical stress/work/heat, coefficient scaling, timestep and
invalid-coefficient checks pass. Additional checks deliberately set molecular
well depths and diameters to zero: prescribed periodic fluxes remain exactly
identical to the original valid-database result, and the prescribed timestep
remains valid and positive. The ordinary molecular path and species-diffusion
path reject the invalid database; negative conserved density also rejects.
The analytical face-flux error remains zero against the unchanged 1e-12 gate.

The first test attempt incorrectly expected a negative temperature *guess* to
invalidate a valid conserved state. Existing temperature inversion legitimately
recovers it. The fixture now uses negative conserved density; no solver guard
was relaxed to satisfy the test. Initial failure is preserved separately.

Provenance, build, syntax, project contract and whitespace checks pass. Existing
real-comparison warnings remain. Pinned YAML SHA-256 remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Artifacts: `build/constant-transport-independent-checks/`, including build,
initial fixture failure, final CTest and `LastTest.log` files. This is focused
kernel evidence, not a new application or physical-validation run.
