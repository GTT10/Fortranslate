# LES SDC supplied halos and coordinated iterations — 2026-09-27

Implementation `74f5c85c`, based on `458523ae`, local main, version 0.246.0.

## Change

The LES SDC wrapper now accepts paired gas/temperature and source-halo callbacks.
Hydro and molecular transport use the existing regular SDC callback path. LES
requests two layers for Smagorinsky/WALE/Vreman, or the Dynamic filter radius
plus five layers. Its ghosted operators consume those supplied values. LES
boundary inventory comes from the exported fluxes on all patch faces, allowing
a future coordinator to distinguish internal interfaces from physical edges.

Optional `new_source_estimate` retains force, LES and optional molecular-transport
families for external iteration scheduling. As with the gas API, the caller
resets gas/temperature to the same interval start before each one-iteration call
and retains the estimate. State, estimate, reaction history, fluxes and diagnostics
publish only after successful completion. The two halo interface types are now
public in the regular reactive module.

## Independently completed checks

Artifacts: `build/les-sdc-halo-checks/` (`provenance.log`, `build.log`, `ctest.log`,
`LastTest.log`). Existing fixed-mechanism Ninja Debug configuration, GNU Fortran
13.3.0; serial focused test.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_sdc_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^les_sdc_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, 11.11 s (11.16 s total), no failures or skips. Provenance, build,
Python syntax, project contract and whitespace passed. Existing warnings remain.

All four closures with molecular transport match the original periodic path
using caller-supplied periodic halos (`1e-12` scaled state criterion). Two
externally coordinated one-iteration calls match a two-iteration reference
from the same zero initial source estimate (`1e-14`). A new-time halo failure
preserves gas, temperature, source estimate, budget and flux outputs. Prior
boundary, inventory and full-face reconstruction checks execute unchanged.

This is a patch-level interface result, not an AMR/MPI integration or reflux
result. Physical-boundary supplied halos, active chemistry, application input,
selected mechanisms, physical validation and external-solver parity were not
newly tested here.
