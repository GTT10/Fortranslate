# Physical EB film level initialization — 2026-09-27

Implementation `7264f48c548c3e300db0f6e878aaac99663893ca`, based on `0eb9d518`,
local main, version 0.246.0.

## Change

`migrate_mpi_amr_wall_films` can now export optional paired `nodes` and `probes`
arrays. They come from the exact ownership calculation used to migrate each film,
including the existing EB surface/normal checks and finest-level selection.
The EB probe is offset into fluid along the film normal; the wall position stays
unchanged. Both metadata arrays publish with the film/component transaction.
Presence must match across ranks, and callers supply both arrays or neither.
The non-EB path reports its existing closed-bound wall-position ownership probe.

`initialize_mpi_eb_film_levels` privately migrates flat input films, then groups
them by their zero-based logical depth into `spray_level_film_population` arrays.
Each group retains films, component masses, node indices and fluid probes in the
same order. Empty levels/ranks have allocated empty arrays. Inputs remain
unchanged, and a failed migration retains any previously supplied hierarchy.

This connects physical film placement to the new SDC population representation.
Virtual/ghost film construction, recursive callback wiring, collision-created
films, outer source synchronization and application routing remain open. It does
not enable recursive films in the application yet.

## Independently completed checks

Artifacts: `build/spray-film-levels-checks/`. Existing fixed-mechanism Ninja Debug
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

**2/2 passed**, 2 ranks 0.53 s and 4 ranks 0.51 s, 1.09 s total.
No build/test failures or skips. Provenance, build, Python syntax, project contract
and whitespace checks passed. Existing compiler warnings remain.

The existing three-level EB film migration fixture now additionally checks:

- every local film appears on exactly one level, with the correct owning rank;
- wall coordinates, normals, materials and component masses are preserved;
- reported probes exactly equal wall position plus the existing fluid offset;
- the covered root and globally empty hierarchy contain correctly shaped empty
  payload arrays;
- grouped global component inventory matches the initial inventory;
- receive-capacity failure retains metadata sentinels and any prior hierarchy;
- rank-disagreed metadata presence rejects collectively.

Existing geometry/normal rejection, owner-map changes, coarsening/refinement
migration and component-inventory checks also passed in these tests. This is
three-level placement evidence, not three-level coupled film evolution. No SDC
gas/film application interval, impact/splash, selected mechanism or physical
validation was run.
