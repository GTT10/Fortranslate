# Recursive MPI EB predictor context — 2026-09-27

Implementation `3008b053`, based on `843d1c52`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`advance_mpi_eb_tree_sdc_subcycle_3d` accepts optional `predictor_builder`
with public interface `mpi_eb_sdc_recursive_predictor_builder`. The adapter
passes combined valid-cell rates, interval, iteration, depth, interval start
relative to the outer step, and the complete owner-local parent predictor
endpoint fields. Parent interpolation fractions equal the midpoint of the
start/end fractions passed to level hydro; nonancestor entries are zero.

The builder returns a complete three-cell predictor halo through the existing
level SDC builder and hydro validation. Callback presence is collectively
agreed. Explicit present/absent call branches preserve the default path and
avoid passing disassociated procedure pointers. The builder works with or
without either auxiliary source callback interface. Its external state must
remain private until the outer step succeeds.

This forwards the context needed for particle source halo assembly. Actual
physical/virtual/ghost population integration, source assembly and recursive
spray application dispatch remain open.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. `spray_mpi_eb_sdc_force_2` passed: 1/1,
36.92 s; CTest total 36.97 s. No build/test failures or skips occurred.

The two-level distributed fixture now builds predictors both with no auxiliary
callback and with the zero-feedback halo-context callback, for densities 2 and
4. It checks six builder calls per outer step, iteration sequence, interval
start/duration, midpoint ancestor fractions and endpoint ownership. The normal
source fill uses the supplied endpoints and fractions. Existing analytic
acceleration/force-density velocity, work and temperature checks pass.

A failure on only the last rank during the sixth builder call (last fine
interval, second iteration) rejects the whole step and preserves gas,
temperature, boundary/wall/chemical/external/auxiliary inventories and limiter.
The existing common-step particle forcing checks also execute.

This focused recursive test disables chemistry/transport and uses zero particle
feedback. Deep hierarchies, actual particle halo modifications and reactive
recursive builders were not independently tested here. No full matrix,
installed audit, external solver numerical parity or physical validation was
performed. Provenance, Python syntax, project contract and whitespace passed.

## Reproduction

Artifacts: `build/sdc-recursive-predictor-checks/`, including build, provenance,
contract, CTest summary and copied detailed `LastTest.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_force_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
