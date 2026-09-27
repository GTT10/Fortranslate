# Periodic film normal probes — 2026-09-27

Implementation `bddb368e`, based on `1cd645b1`, local main, version 0.246.0.

## Change

Film surface validation now splits a short normal probe at periodic faces and
passes each root-contained segment to the existing tetrahedral intersection.
The returned contact is lifted back to the original wall coordinate frame.
Ownership permits exterior probe coordinates only on periodic axes and wraps
the cell lookup. Film payload coordinates still have to lie inside the closed
physical domain and remain unchanged. Nonperiodic validation uses the original
intersection directly. Probes longer than a domain length reject; the split
loop is bounded to four segments. Existing collective geometry, capacity and
transactional publication checks remain in place.

## Independently completed checks

Artifacts: `build/spray-periodic-film-probe-checks/` (`provenance.log`,
`build.log`, `ctest.log`, `LastTest.log`). Existing fixed-mechanism Ninja Debug
build with GNU Fortran 13.3.0 and OpenMPI 4.1.6, `/usr/bin/mpifort` and
`/usr/bin/mpiexec`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_eb_film test_spray_mpi_eb_reactive --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_eb_(film_(2|4)|sdc_recursive_impact-periodic_2)$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**3/3 passed**, zero failures or skips: film migration/copies on two ranks
0.49 s, four ranks 0.51 s, periodic recursive impact coupling 8.39 s; total
9.44 s. Provenance, build, Python syntax, project contract and whitespace checks
passed. Existing compiler warnings remain.

The new fixture uses a periodic piecewise-linear x field plus an oblique y
surface, with matching slope at the x seam. Films at x=0 and x=4 exercise an
exterior solid probe and exterior fluid probe respectively. It checks preserved
wall endpoints, unwrapped fluid probes, assignment to the left fine patch,
physical hierarchy construction and both root virtual copies. Existing film
rollback, capacity, component inventory and nonperiodic cases run unchanged.

The recursive coupled test is an existing planar-wall regression: it does not
establish thermal evolution of the new oblique seam fixture. No new application,
restart, selected-mechanism, physical-validation or external-parity run was made.
