# Periodic SDC reaction coupling — 2026-09-26

Implementation `7a94c85d`, parent `5cbfac80`, local main, version 0.246.0.

The periodic SDC library API accepts an optional mechanism and persistent
species/energy reaction-source array. Previous reaction rates participate in
hydro prediction and the state used to evaluate new external sources. Each
iteration then rebuilds nonreacting rates, integrates each cell from its
original species densities for the full interval, and updates reaction rates
after subtracting the nonreacting contribution. Internal-energy forcing is
derived from total energy minus kinetic energy at the nonreacting endpoint.
Density follows the integrated species sum; momentum retains its nonreacting
update. Source, state and temperature outputs are published only on success.

The ordering follows `Advance.cpp::do_sdc_iteration`, `Hydro.cpp` and
`React.cpp::react_state` from pinned PeleC
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`, inspected with `git show` in
`/mnt/sn850/Projects/gtt2/external/PeleC`. Integration uses PeleF's adaptive
implicit reactor, not a claim of identical reference reactor numerics.
Callers can retain the optional reaction-source array between steps; absent
that array, each call starts its reaction-source estimate at zero.

## Completed checks

- GNU Fortran 13.3.0 Debug, fixed tests-enabled build
  `build/spray-amr-reactive-debug`: target `test_reactive_ppm_step_3d` passed.
- `unit_reactive_ppm_step_3d`: **1/1 passed**, 19.96 s test / 20.04 s total.
  Added a uniform 4³ reacting case with two SDC iterations, compared against
  one full-interval cell integration. Checks require an active species change,
  matching species/temperature, reaction-source accounting, and unchanged
  state/temperature/reaction source after an external-source callback failure.
  Existing nonreacting SDC and hydro/source tests remain active.
- Initial fixture had no radical seeds and failed the active-reaction gate.
  It was replaced by a mixture with H/O/OH seeds; no tolerances were relaxed.
  A CTest launched before build completion is retained separately and is not
  counted. Final evidence is the post-build run above. No final failures/skips.
- Provenance, Python syntax, project contract and whitespace passed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_step_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_step_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/sdc-reaction-checks/` contains build, provenance, CTest,
LastTest and initial fixture logs. Fixed/selected application input still
rejects reacting SDC; application dispatch, nonuniform reacting qualification,
physical boundaries, AMR/MPI and external solver comparison remain open.
