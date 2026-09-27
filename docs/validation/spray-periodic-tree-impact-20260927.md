# Periodic atomic EB tree impacts — 2026-09-27

Implementation `500ab29e`, based on `7baf16a5`, local main, version 0.246.0.

## Change

The standalone `advance_mpi_eb_tree_impacts` API now accepts optional
`boundary_faces`, consistent with its existing level-set flight kernel. Explicit
faces require the outflow ledger and valid paired boundary labels; inflow is
rejected. Periodic axes reach initial/final film migration and final particle
migration, and all boundary labels reach flight. Omission preserves the existing
outflow behavior. Pressures stay attached to original parcels until flight;
state, films, IDs and ledgers publish only after every migration succeeds.

## Independently completed checks

Artifacts: `build/spray-periodic-tree-impact-checks/` (`provenance.log`,
`build.log`, `ctest.log`, `LastTest.log`). Existing fixed-mechanism Ninja Debug
build, GNU Fortran 13.3.0, OpenMPI 4.1.6, `/usr/bin/mpifort`, `/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_impact --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_impact_(2|4)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, 0.50 s each, 1.05 s total; no failures or skips. Provenance,
build, Python syntax, contract and whitespace passed. Existing warnings remain.

The added fixture crosses the periodic z face before depositing on an oblique
x/y wall. It checks wrapped contact z=0.09, finest film ownership, zero exterior
loss and combined wall/particle/film inventory at the existing `3e-12` criterion.
A late film receive-capacity failure preserves original positions, populations,
ID watermark, hit count and both ledgers; retry succeeds. Existing splash,
outflow, migration and rejection tests execute in the same programs.

This is a focused library result. No new application, thermal-evolution,
selected-mechanism, physical-validation or external-parity run was made.
