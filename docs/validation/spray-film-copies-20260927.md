# Virtual and ghost EB film populations — 2026-09-27

Implementation `f87992719f97849fad643f4631a1f97d5fe78649`, based on `a4b9edfc`,
local main, version 0.246.0.

## Change

`refresh_mpi_eb_film_copies` builds two private outputs from the physical film
hierarchy without modifying it:

- virtual films represent every descendant physical level, migrated to the
  requested parent level using the existing capped EB ownership calculation;
- incoming ghost films represent only the current level's physical films that
  lie in the next level's valid or grown boxes, located using their fluid probes.

The ghost lookup shares the current particle convention: valid boxes take
precedence, then grown boxes, with tree-order ties and at most one destination
per input film. It does not add periodic image routing or virtual aggregation.
Wall coordinates, normals, temperature, geometry, material, component masses and
fluid probes remain exact payloads. Receive order is source rank/local order.

Physical level/owner associations and probe consistency are validated. MPI ranks
agree on copy controls, and existing migration validates replicated topology and
EB geometry before ghost communication. Both virtual and ghost outputs publish
only after both phases succeed. A late ghost-capacity failure therefore retains
the previous virtual output as well. Root-empty and finest-level outgoing sets
use allocated empty arrays.

This completes film copy construction as a library operation. Recursive clock
coordination, gas/source callback integration, collision-created films, outer
publication and application routing remain required. No recursive film mode is
enabled by this change alone.

## Independently completed checks

Artifacts: `build/spray-film-copies-checks/`. Existing fixed-mechanism Ninja Debug
build `build/spray-amr-reactive-debug`; GNU Fortran 13.3.0 and OpenMPI 4.1.6
(`/usr/bin/mpifort`, `/usr/bin/mpiexec`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_film --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_film_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**, 0.50 s on each rank layout, 1.06 s total. No build/test failures
or skips. Provenance, build, Python syntax, contract and whitespace checks passed.
Existing compiler warnings remain.

The three-level EB fixture verifies root virtuals from all descendants,
intermediate virtuals from the nested level, and nearby parent films delivered
to the child's grown box. It checks exact payloads, target nodes and probes,
unchanged physical populations and empty finest-level outgoing copies. A larger
grown box makes the ghost receive exceed capacity after the virtual receive has
succeeded; both prior outputs remain unchanged. Rank-disagreed growth controls
also reject collectively. The existing migration and level-placement checks pass
in the same two tests.

These are three-level copy-placement checks, not coupled film evolution.
No SDC gas/film interval, impact/splash case, application runtime, selected
mechanism or physical validation was run.
