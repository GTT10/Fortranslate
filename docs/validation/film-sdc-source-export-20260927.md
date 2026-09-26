# MPI AMR/EB film source export — 2026-09-27

Implementation `7e2e2a97`, based on `1bfdabe4` on local `main`.
PeleC reference `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Change

`advance_mpi_amr_film_coupling` accepts optional `gas_change`, with the same
owner-local, time-integrated conserved source semantics as the particle interval.
Each successful gas deposition exports its independent contribution; accepted
substeps add it to a private zero field with overflow checks. Failed trials,
substep exhaustion and final film migration failures do not publish the source.
Successful publication replaces the previous field, alongside gas, temperature,
film state, component masses and wall heat. Output presence is agreed before
collective operations. Hidden parents follow the existing deposition restriction;
nonowners retain unallocated arrays. Empty populations and disabled evolution
produce zero source. The gas background is never subtracted to infer feedback.
Wall heat remains in the separate boundary ledger.

This is a prerequisite for film SDC scheduling. It does not enable film impacts
or film evolution in the SDC wrapper/application; their existing rejection remains.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed H2O2 tests-enabled build
`build/spray-amr-reactive-debug`; `/usr/bin/mpifort` and `/usr/bin/mpiexec`,
OpenMPI 4.1.6. Builds completed; existing compiler warnings remain.
CTest **2/2 passed**, 1.24 s total:

- `spray_mpi_eb_film_step_2`: 0.63 s.
- `spray_mpi_amr_film_step_2`: 0.56 s.

Both use three-level trees, two ranks, mixed water/peroxide films, evaporation,
heating and wall energy. Exported field ownership matches gas ownership and
reconstructs gas changes within 3e-13 times max(1,abs(final component)).
The existing combined gas/film/wall conservation bound of 2e-12 times
max(1e-12,abs(initial inventory)) passes. Substep exhaustion preserves sentinel
source values as well as film and gas state. The EB test additionally rejects
rank-disagreeing source/geometry presence and nonuniform wall temperature.
Globally empty films replace the previous source with zero.

Provenance, Python syntax, project contract and diff whitespace passed.
No build/test failures occurred. These focused checks do not establish a full
configuration matrix, selected-mechanism qualification, installed-binary audit,
physical validation, order accuracy or upstream solver parity.

## Reproduction

Artifacts: `build/film-sdc-source-checks/` (`provenance.log`, `build.log`,
`ctest.log`, `detail.log`, `contract.log`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_film_step test_spray_mpi_eb_film_step --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_(amr|eb)_film_step_2$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
