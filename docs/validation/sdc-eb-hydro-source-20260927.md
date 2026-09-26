# Conserved SDC sources in EB CTU prediction — 2026-09-27

Implementation `f607bb1e`, base `6d386978`; PeleF 0.246.0.
PeleC pin `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

## Implementation

`compute_reactive_eb_plm_ctu_fluxes_3d` accepts optional `conserved_source`,
mutually exclusive with primitive source/energy inputs. The existing
`reactive_pelec_source_to_primitive` conversion produces primitive and internal
energy rates at active cells; covered cells are skipped. The existing final
two-direction transverse correction adds these rates with dt/2. Normal and
intermediate corner flux stages remain source-free. Candidate predictions are
published only after all conversion and flux stages succeed.

Pinned `Source/Godunov.cpp` lines 908–913 were inspected with `git show`: the
EB final X correction calls `pc_transdd` with `srcQ` and `hdt` for uncovered
cells. The existing Fortran transverse routine supplies that stage ordering;
this change connects the conserved-rate input. No new compiled C++ oracle was
run, and this is not a full EB SDC step or full solver parity claim.

## Focused evidence

Artifacts: `build/sdc-eb-hydro-source-checks/` (local). GNU Fortran 13.3.0,
existing Ninja Debug tests-enabled build `build/spray-amr-reactive-debug`;
compiler/backend unchanged.

- Provenance passed: `python3 tools/check_mechanism_provenance.py --bundle
  mechanisms/h2o2_full.json`; SHA
  `0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
- Build passed: `env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build
  build/spray-amr-reactive-debug --target test_reactive_eb_hydro_3d
  --parallel 2` (`build.log`). Compiler warnings remain.
- CTest **1/1 passed**, 25.36 s (25.41 s total): `env -u LD_LIBRARY_PATH
  -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R
  '^regression_reactive_eb_hydro_3d$' --output-on-failure` (`ctest.log`).
- Extended 8-cubed regular and x/y/z plane-cut working geometries. A conserved
  momentum/work source corresponding to x acceleration 8 matches explicitly
  supplied primitive acceleration fluxes within 2e-12 relative scale in all
  directions; the nonzero source changes the source-free flux. Covered source
  payloads contain NaNs and are ignored. Invalid active source rejects before
  publication; existing failed-stage rollback checks retain predictions.
- Existing normal/transverse, CTU, physical boundary and redistribution checks
  in the regression executable also passed.
- Python syntax for `tools/check_project_contract.py`, project contract and
  `git diff --check` passed.

No coupled EB SDC, application run, MPI path, selected build, full matrix,
installed audit or physical validation was run. The full port remains incomplete.
