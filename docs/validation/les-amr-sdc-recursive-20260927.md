# Recursive AMR LES SDC — 2026-09-27

Implementation `0313121d`, based on `65db541a`, local main, version 0.246.0.

## Change

The AMR SDC diffusion callback now receives coarse old/new packed endpoints in
addition to level selection and stage fractions. The recursive driver forwards
it to every level interval. LES RHS evaluation selects that depth, retains zero
rates elsewhere, replaces lower-level donors by interval-start endpoints and
passes future endpoints/fractions into existing temporal coarse ghost filling.
Endpoint shapes, presence and fraction bounds are validated before evaluation.

`advance_amr_patch_tree_les_sdc_3d` now accepts `subcycling=.true.` and optional
`level_steps`. Common-step behavior remains the default. Both modes use the same
LES stage evaluator; the existing recursive SDC driver integrates local-interval
flux density, accumulates child fluxes and synchronizes with the parent.
LES also participates in parent old/new predictor source fields. Caller state,
budgets, diffusivity and step counts are published only on success.

## Independently completed checks

Artifacts: `build/les-amr-sdc-recursive-checks/` (`provenance.log`, `build.log`,
`ctest.log`, `LastTest.log`). Existing fixed-mechanism Ninja Debug configuration,
GNU Fortran 13.3.0. Both selected tests are serial.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_amr_sdc --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^les_amr_sdc(_recursive)?$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**2/2 passed**: common-step 10.31 s, recursive 15.36 s, total 25.72 s. No failures
or skips. Provenance, build, Python syntax, contract and whitespace passed.
Existing compiler/internal-callback warnings remain.

Both modes run Smagorinsky, WALE, Vreman and Dynamic on root-only and two-level
periodic trees. Root-only states match the regular LES SDC reference within
`1e-12`. Recursive counts are exactly one root interval and two fine intervals;
common-step counts are one each. Three models include molecular transport, while
Vreman exercises the LES-only diffusive path. Composite inventory residual is
zero at the printed precision for this fixture, satisfying `3e-12` after reflux.
Invalid filter input preserves all patch states, temperatures, budgets,
diffusivity and step-count outputs.

This does not establish temporal convergence, active chemistry, physical-boundary
AMR, deeper hierarchies, MPI, application input, selected-mechanism coverage,
physical validation or external-solver numerical parity.
