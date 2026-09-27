# Common-step AMR LES SDC — 2026-09-27

Implementation `ef7bd8c4`, based on `739ddbdd`, local main, version 0.246.0.

## Change

The AMR SDC driver accepts an optional diffusion evaluator returning patch rates
and instantaneous face flux density. It validates returned shapes/finiteness,
adds rates to the diffusive source used by predictors and corrections, multiplies
face fluxes by dt, and combines old/new fluxes with hydro and molecular transport
before existing reflux and restriction. The callback also works when molecular
transport is disabled. Existing private-state publication preserves rollback.

`compute_amr_patch_tree_les_rhs_3d` evaluates all patches from one immutable stage
state, using existing coarse/fine ghost filling and four LES closures. It returns
unrefluxed rates and flux density. `les_amr_sdc_mod:advance_amr_patch_tree_les_sdc_3d`
connects this evaluator to common-step SDC, including existing gas chemistry,
transport and physical-boundary arguments. LES uses diffusive flux accounting;
it is not included in the volume-force ledger.

This wrapper advances every level with a common interval. It does not yet expose
recursive temporal donors, MPI ownership or application input selection.

## Independently completed checks

Artifacts: `build/les-amr-sdc-checks/` (`provenance.log`, `build.log`, `ctest.log`,
`LastTest.log`). Existing fixed-mechanism Ninja Debug build, GNU Fortran 13.3.0.
The focused test is serial and links the opt-in tree runtime.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_amr_sdc --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^les_amr_sdc$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, 10.22 s (10.27 s total), no failures or skips. Provenance, build,
Python syntax, project contract and whitespace passed. Existing warnings and
internal-callback executable-stack warnings remain.

The test exercises Smagorinsky, WALE, Vreman and Dynamic on a periodic root-only
tree and a two-level tree with a half-domain fine patch. Root-only gas states
match the regular LES SDC API within `1e-12`. Three models include molecular
transport; Vreman exercises LES without molecular transport. The normalized
composite inventory residual after combined reflux is zero at the printed
precision for this fixture, satisfying the `3e-12` criterion. An invalid Dynamic
filter preserves every patch state/temperature and the supplied budget/diffusivity
outputs.

No active chemistry, physical-boundary AMR, deeper hierarchy, recursive SDC, MPI,
application, selected-mechanism, physical-validation or external-solver parity
run was performed here.
