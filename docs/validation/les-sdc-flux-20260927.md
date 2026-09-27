# Combined LES SDC face flux export — 2026-09-27

Implementation `d01fbcf3`, based on `821e72af`, local main, version 0.246.0.

## Change

Regular algebraic LES and Dynamic boundary RHS APIs now optionally export all
three full-face flux arrays, with shape validation and publication after success.
The periodic algebraic kernel uses physical cell indices and duplicates its
last periodic face at index zero; ghosted physical-boundary exports retain
their active-cell offset convention.

`advance_les_sdc_3d` now accepts optional `flux_x`, `flux_y`, `flux_z` together.
It retains the old LES flux and the final new-stage flux, adds their half-weighted
average to the SDC hydro/transport flux, and publishes all outputs with the
accepted state. Arrays have one extra face in their normal direction. Values
are time-averaged flux density in the positive coordinate direction: multiply
by interval and face area for an integrated flux. This supplies the data for
AMR synchronization; it does not yet connect an AMR LES SDC driver.

## Independently completed checks

Artifacts: `build/les-sdc-flux-checks/` (`provenance.log`, `build.log`, `ctest.log`,
`LastTest.log`). Existing fixed-mechanism Ninja Debug configuration, GNU Fortran
13.3.0; the focused test is serial.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_sdc_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^les_sdc_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, 7.35 s (7.40 s total), no failures or skips. Provenance, build,
Python syntax, project contract and whitespace passed. Existing compiler and
internal-callback executable-stack warnings remain.

For all four closures, periodic and isothermal-wall cases with molecular
transport reconstruct every updated conserved cell from the initial state
minus the divergence of exported combined fluxes (`1e-12` scaled criterion).
These reconstruction cases omit volume forcing/reactions so that the full
update is a flux divergence. The preceding nonzero-force inventory and zero-LES
checks still execute. An invalid Dynamic filter and malformed face-array shape
preserve state and caller flux outputs.

No AMR reflux, MPI, application, active chemistry, selected-mechanism matrix,
physical-validation or external-solver parity run was performed here.
