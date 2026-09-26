# Recursive MPI EB gas acceleration — 2026-09-27

Implementation `66061f0c`, based on `dba81a24`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Change

`advance_mpi_eb_tree_sdc_subcycle_3d` now accepts optional
`gas_acceleration(3)` and forwards it to each selected-level SDC interval,
with or without an auxiliary source callback. Omission retains zero acceleration.
The existing level driver validates finite, bounded, rank-consistent values and
forms the momentum source as `force + density*gas_acceleration`, with work from
the local gas velocity. Recursive candidate state and ledgers remain private
until the outer step succeeds.

This connects the gas forcing prerequisite for recursive spray coupling.
Recursive physical/virtual/ghost particle scheduling, film evolution and
application dispatch remain unconnected. The serial recursive gas API is not
extended by this change.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. The existing `spray_mpi_eb_sdc_force_2`
test passed (1/1, 24.62 s; CTest total 24.66 s).

The test now runs its two-level, non-unit-density gas fixture through both the
common-step spray wrapper with empty particles and the recursive gas driver.
For densities 2 and 4, acceleration [1000,0,3000], force density [20,0,40],
initial velocity [1,0,3] and interval 1e-8, it checks velocity increment
`dt*(acceleration+force/density)` within 1e-12, the integrated impulse/work
within the existing relative tolerance 3e-13, and unchanged temperature within
1e-8. The new recursive calls use no auxiliary source callback. Existing
particle acceleration assertions execute in the same test.

Provenance preflight, target build, Python syntax, project contract and
`git diff --check` passed. No build/test failures or skipped tests occurred.
This focused test does not establish a full configuration matrix, installed
binary qualification, recursive particle coupling, temporal order, or physical
validation.

## Reproduction

Artifacts: `build/sdc-recursive-acceleration-checks/` (provenance, build, tests,
project contract and copied CTest detail logs).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sdc_force_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
