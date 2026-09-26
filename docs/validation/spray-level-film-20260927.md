# Liquid-film coupling in capped AMR views — 2026-09-27

Implementation `5f1d5a21`, based on `2e9d3cf2`, local `main`, version 0.246.0.
Pinned PeleC reference: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.
H2O2 provenance SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

## Change

Optional `maximum_level` caps film ownership in `migrate_mpi_amr_wall_films`.
MPI configuration agreement includes the effective ceiling; the original tree
and EB nodal-geometry checks remain in place. Exact wall coordinates and all
film/component fields retain the existing transactional exchange semantics.

The film thermal step forwards this ceiling to initial/final migration, gas
sampling and deposition, and uses capped spacing for planar wall gas probes.
The particle step forwards the ceiling to impact-film migration and therefore
no longer rejects capped calls merely because film populations are supplied.
An explicit ceiling returns a provisional gas view without parent restriction;
the caller must synchronize sources and keep physical/virtual inventories
separate. Recursive scheduling and application integration remain open.

## Completed verification

Fixed-mechanism Ninja Debug build `build/spray-amr-reactive-debug`, GNU Fortran
13.3.0, OpenMPI/OpenRTE 4.1.6. Four tests passed in 2.86 s total:

- `spray_mpi_amr_film_step_2`: 0.59 s.
- `spray_mpi_amr_step_2`: 0.87 s.
- `spray_mpi_amr_film_step_4`: 0.57 s.
- `spray_mpi_amr_step_4`: 0.78 s.

The new film fixture caps a three-level tree at level zero and compares it
with a root-only tree. Both contain films at the two exact physical x walls.
Film ownership/count, components, positions, normals, mass, temperature,
diameter, height, material IDs, root gas/temperature, substep count and wall
inventory match exactly. Finer gas and temperatures remain unchanged and
their feedback is zero. Rank-disagreeing ceilings reject, preserving gas,
components/count and the caller's wall ledger. Existing film conservation,
heating/evaporation and late rollback assertions also execute, together with
the particle coupling regression tests.

The newly capped EB-film and impact-film paths were compiled but not separately
exercised. No full matrix, installed audit, recursive application, external
solver parity or physical validation was performed. Provenance, target build,
Python syntax, project contract and whitespace checks passed. There were no
build/test failures or test skips in this change.

## Reproduction

Artifacts: `build/spray-level-film-checks/`, including copied CTest details.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_film_step test_spray_mpi_amr_step --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_(film_step|step)_[24]$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
