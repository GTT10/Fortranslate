# Recursive SDC source access to ancestor halo context — 2026-09-27

Implementation `46045464`, based on `796048be`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Interface

The MPI recursive EB SDC driver accepts `halo_source`, with interface
`mpi_eb_sdc_recursive_halo_source`, as an alternative to `auxiliary_source`.
Both use the existing `auxiliary_change` inventory and optional cross-level
source routing. Simultaneously supplying both callbacks rejects. MPI agreement
distinguishes the callback signatures before the recursive loop.

The new callback receives the existing state, temperature, interval, stage,
depth and interval-start time, plus read-only ancestor state/temperature
endpoints and interpolation fractions. Stage zero selects interval-start
fractions; positive stages select interval-end fractions, using the recursive
ancestor clocks. Nonancestor entries are zero. These inputs can feed the
stage-time particle gas halo operator without reconstructing clocks separately.
Callback-owned particle state must remain private until outer acceptance.
This interface does not yet implement recursive population scheduling.

## Completed verification

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. `spray_mpi_eb_sdc_force_2` passed
(1/1, 33.44 s; CTest total 33.49 s).

The existing two-level accelerated gas fixture now runs a third mode through
the new callback at densities 2 and 4. The callback returns zero source and
checks the fraction vector shape, old/new fractions against interval timing
within 1e-14, owner-local snapshot allocation and state/temperature shapes.
During fine steps it verifies that the ancestor endpoint contains the completed
coarse momentum update. Existing analytic velocity, forcing impulse/work and
temperature assertions still pass. Common-step and recursive calls without
this callback remain exercised by the other two modes.

Provenance, target build, Python syntax, project contract and whitespace checks
passed. No build/test failures or skips occurred. The new callback was checked
with zero feedback on two levels, not with particle evolution, cross-level
feedback, deeper hierarchies or chemistry. No full matrix, installed audit,
recursive spray application, temporal-order study, external solver parity or
physical validation was performed.

## Reproduction

Artifacts: `build/sdc-particle-halo-context-checks/`, including CTest details.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_force_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
