# EB flux and Euler viscosity dispatch — 2026-09-25

## Source and scope

Code `2e83d78`, parent `500c255`. The EB Riemann-flux API accepts optional
`difmag`, a volume-fraction halo and physical no-slip masks. Positive
coefficients apply the masked viscosity kernel to candidate intensive fluxes
before publication. Face area/aperture multiplication remains in the existing
geometry divergence/integrated-flux path.

Without an external halo, the routine extends state and volume fraction by
outflow extrapolation. With an external state halo, positive coefficients
require a matching one-cell volume-fraction halo, including edges/corners;
its interior must equal the supplied geometry. Interior states come from the
primary state array. Zero/omitted coefficients retain the old face-only halo
contract and zero correction policy. Invalid input publishes no flux.

The conservative RHS and plain Euler API forward all three arguments. Flux
redistribution and state redistribution Euler APIs forward an optional
coefficient through their existing outflow path. Existing state and integrated
flux transaction boundaries remain. EB tree/MPI halo producers, higher coupled
drivers, input/default/restart policy remain outstanding; this is not a full
EB solver port.

## Completed checks and build failure

GNU Fortran 13.3.0, Ninja Debug, Python 3.12.3; reused
`build/spray-amr-reactive-debug`. Focused CTest **1/1 passed**, 0.95 s
(0.99 s elapsed). Initial build failed because a new argument line exceeded
132 columns; wrapping it fixed the build. No runtime failures or skips.
Compiler warnings remain.

The axis-plane cut-cell fixture now uses density/velocity variation and
coefficient 0.037. It checks changed open flux, zero closed flux, plain Euler
agreement with divergence of corrected intensive flux (scaled 2e-14), success
of both redistribution paths, and negative-coefficient state/temperature
rollback. Covered states are zero in the fixture. Existing uniform tangential,
small-cell and conservation checks also pass. External positive-coefficient
halos and no-slip masks were not separately exercised by this test. No
tolerances were loosened. This is focused API evidence, not external PeleC
parity or physical validation.

Mechanism provenance, Python syntax, project contract and whitespace checks
passed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_eb_hydro_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^regression_reactive_eb_hydro_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-viscosity-dispatch-checks/` contains failed `build.log`,
successful `rebuild.log` and copied `LastTest.log`.
