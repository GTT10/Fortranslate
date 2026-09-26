# Particle predictor halo composition — 2026-09-27

Implementation `ddecc1e1`, based on `32183e2f`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Implementation

`fill_sdc_mpi_eb_tree_source_halo_3d` adds optional `same_level_only`. This
mode fills valid/peer/periodic donors while leaving physical exterior and
missing same-level donors zero. MPI agrees the mode, and ancestor interpolation
arguments are rejected in this mode. The default fill remains unchanged.

`fill_sdc_mpi_eb_particle_predictor_halo_3d` composes a complete signed predictor
from total valid rates, their particle component and the unmatched deposited
particle halo remainder. Ghost cells use:

```text
normal total halo - normal particle halo + peer particle halo + particle remainder
```

This removes particle contributions introduced by normal parent/boundary
extrapolation before adding the actual deposited halo. Nonparticle contributions
retain normal interpolation. Total/particle ancestor endpoints and their shared
fractions are supplied together or omitted together. The caller must supply the
particle component actually included in the total rates and a same-stage,
same-unit unmatched remainder from the source transfer.

Valid-cell total rates are retained exactly. The remainder must have owner-local
halo shape, finite values and zero valid-cell values; nonzero covered ghost
remainder is also rejected. Spatial lower-bound offsets are accepted. Signed
arithmetic is overflow-checked; collective failures preserve the prior output.
This adds halo assembly, not the full recursive spray wrapper or application
connection. It does not validate the caller's source-family provenance.

## Completed checks

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. `unit_mpi_eb_tree_ctu_hydro_3d_3` passed:
1/1, 25.74 s; CTest total 25.79 s. No build/test failures or skips occurred.
Provenance, Python syntax, project contract and whitespace checks passed.

The three-level EB fixture runs distributed and all-on-rank-zero layouts. Its
particle-only predictor verifies that the composed result equals peer particle
sources plus a supplied unmatched remainder, with no duplicated ancestor or
physical-boundary extrapolation. Temporal endpoint arguments are exercised with
a shared 0.5 fraction. The composed halo is then accepted by actual SDC hydro.
Adding a remainder inside valid cells rejects and retains every previous halo
array. Existing temporal-source serial/MPI parity, source-aware hydro, supplied
halo sensitivity and rollback checks also execute.

The composition fixture uses identical total/particle source families and
constant temporal endpoints. Mixed source-family temporal differences and
physical-layout remainder bounds were not independently exercised here. No
recursive gas/particle run, full matrix, installed audit, external solver
numerical parity or physical validation was performed.

## Reproduction

Artifacts: `build/spray-predictor-composition-checks/`, including initial/final
build, provenance, contract, CTest summary and copied detailed `LastTest.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_eb_tree_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_mpi_eb_tree_ctu_hydro_3d_3$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
