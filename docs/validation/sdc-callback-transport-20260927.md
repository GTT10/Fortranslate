# SDC additional sources with molecular transport — 2026-09-27

Implementation `43c71c21`, based on `ef060d4e`, local main, version 0.246.0.

## Change

The regular 3D SDC driver no longer rejects a source evaluator when molecular
transport is present. `source_count` includes a final transport family. The
callback receives only the preceding families, then the driver computes and
stores molecular flux divergence in the final slot. At least one callback
family is required when combining the two. Without transport, the callback
receives the same full source array as before. The builtin force/transport path
retains its existing one/two-family layout.

Old/new source assembly now uses the same dispatcher in both modes. Transport
limiting, physical-boundary inventories, full-face flux export, reaction-source
handling and optional external iteration estimates retain their existing paths.
This removes a prerequisite blocker for LES/additional-source integration; it
does not yet connect an LES closure or add a new application input.

## Independently completed checks

Artifacts: `build/sdc-callback-transport-checks/` (`provenance.log`, `build.log`,
`ctest.log`, `LastTest.log`). Existing fixed-mechanism Ninja Debug configuration,
GNU Fortran 13.3.0. The selected test is serial and links `pelef_core`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_ppm_step_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_ppm_step_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, 23.34 s (23.39 s total), no failures or skips. Provenance, build,
Python syntax, project contract and whitespace passed. Existing warnings remain.

The new test splits a known force into two callback families (one quarter and
three quarters) with molecular transport in a third family. Periodic and x-outflow
cases compare against the existing builtin force-plus-transport path: state and
temperature within `3e-14`, source/boundary inventories and exported face fluxes
within `1e-12` scaled criteria. It verifies one old and two new callback calls,
and checks exact state/temperature/iteration-estimate rollback after a new-time
callback rejection. Existing PPM, SDC source iteration, reactions, boundaries,
transport, halos and flux reconstruction checks execute in the same program.

No LES integration, new application run, MPI/AMR matrix, physical validation or
external solver comparison was performed in this increment.
