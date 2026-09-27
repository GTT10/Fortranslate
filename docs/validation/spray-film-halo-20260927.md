# Frozen stage-time EB film exchange — 2026-09-27

Implementation `ed1143cf925d89156334614c8a94cbbca8092fd5`, based on `2602ee1e`,
local main, version 0.246.0.

## Implemented connection

`kick_mpi_eb_film_halo` supplies the film exchange kernel needed by recursive
SDC. It reads caller-provided stage gas and temperature halos at fluid-side
probes, advances the existing mixed-film frozen heat/evaporation model for the
supplied interval, and returns time-integrated gas exchange and wall energy.
The stationary wall position is never replaced by the fluid probe position.

Deposition uses the containing-cell rule with cut-cell fluid volume. The result
is split into same-level valid source and retained source halo; no copy source
is automatically sent to a coarse level. Optional `physical_source` deposits the
actual exchange at the finest containing composite leaf for later conservative
synchronization. Callers must request this export only for physical films.

Film payloads, component masses, source, remainder, optional physical export
and wall ledger publish together after collective success. MPI controls include
interval, level, halo width, component count, physical-export presence and wall
heat/evolution controls. Model tables and material options remain identical
across ranks by caller contract, as in the existing coupling APIs.

This is a library connection, not completed recursive application support.
Level/copy film populations, old/new stage caching, collision-created films,
outer publication and application routing still need integration. The kernel
expects valid live input films; callers remove depleted outputs before another
call. It does not internally subcycle or change gas state: an inadmissible half
step rejects so its enclosing transaction can retry. Existing stationary-film
models and empirical correlations are unchanged.

## Independently completed checks

Artifacts: `build/spray-film-halo-checks/`. Existing fixed-mechanism Ninja Debug
build `build/spray-amr-reactive-debug`; GNU Fortran 13.3.0, `/usr/bin/mpifort` and
`/usr/bin/mpiexec` OpenMPI 4.1.6.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_sample --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_sample_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 tests passed**, 0.78 s on 2 ranks and 0.79 s on 4 ranks, 1.61 s total.
No test failures or skips. Provenance, final build, Python syntax, contract and
whitespace checks passed. The initial build used an incorrect keyword for the
source-transfer routine; it was corrected to `halo_volume_fraction` before tests
ran. The failed build is preserved in `build-keyword-failed.log`.

The added planar EB fixture has two stationary films on the fine logical level,
one inside its valid cells and one outside in its halo. It verifies:

- component loss, temperature and wall heat against the existing local frozen
  film model using the same uniform gas state;
- exact stationary wall coordinates;
- integrated optional physical exchange on finest composite leaves against the
  extensive film exchange (including the film outside the logical fine patch);
- film internal energy plus gas exchange equals heat supplied by the wall;
- after material 1 updates privately, material 2's temperature falls outside its
  configured conductivity range: every film payload, source field, remainder,
  physical export and wall ledger retains its prior value.

No coupled recursive-film interval, impact/splash case, depletion continuation,
application runtime, selected mechanism or physical validation was run.
