# MPI AMR gas/particle/film remap — 2026-09-18

Implementation: `3ddfa2e8cb64cca8389d2b84e6991867b570ada0`, based on `90444a1`.

`remap_mpi_amr_spray` accepts optional stationary film arrays, component masses
and receive capacity as a group. Film presence agrees collectively before the
optional collective path. Private film migration follows private particle and
gas remapping. Only after all succeed are the topology, owner map, gas,
temperatures, particles/components and films/components published together.
No-film calls keep the existing remap path.

## Completed focused verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja, `/usr/bin/gfortran`,
`/usr/bin/mpifort`, `/usr/bin/mpiexec`; reused unchanged fixed tests-enabled and
selected H2/O2 tests-disabled builds.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_remap --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_remap_(2|4)$' --output-on-failure --output-log build/spray-amr-film-remap-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, zero failures/skips: two ranks **0.52 s**, four ranks **0.55 s**,
total **1.12 s**. The existing three-level nonuniform gas/particle fixture now
includes an upper-wall mixed film. Coarsening, moving refinement to the upper
wall and returning to the original topology preserve all film fields and
component masses exactly and place the film on the expected finest owner.
Existing serial gas remap equality, composite gas inventory conservation and
particle payload checks remain active.

A zero film capacity rejects after private gas/particle remapping and preserves
original topology, owners, gas, temperature, particles and films. Existing
particle-capacity, invalid-source-gas and incompatible-root rollback checks
also run with film state present.

Fixed and selected builds, provenance, syntax, contract and whitespace checks
passed. No build/test attempt failed; existing compiler warnings remain. Selected
evidence is compilation only. No Python source changed.

Artifacts: `build/spray-amr-film-remap-checks/{build,build-selected,ctest}.log`.
Executable: `build/spray-amr-reactive-debug/tests/test_spray_mpi_amr_remap`.

This connects supplied films to the library remap. Impact creation, film-aware
load balancing, application/restart film state, EB interaction and physical
validation remain open. No full matrix, convergence or scaling claim is made.
