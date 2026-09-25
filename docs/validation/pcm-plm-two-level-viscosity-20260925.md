# Serial two-level PCM/PLM viscosity — 2026-09-25

## Source and scope

Code `00d6a7d`, parent `bd40713`. Serial two-level AMR accepts positive
PCM/PLM coefficients. Periodic coarse stepping with an explicit positive
coefficient uses the regular boundary hydro dispatcher; physical coarse
stepping already forwards it. Each fine SSPRK2 stage builds the corresponding
coarse-time-interpolated ghost state and applies the existing regular-cell
artificial-viscosity correction before divergence and mean flux publication.
Reflux therefore includes the corrected coarse and fine fluxes. Physical wall
masks are intersected with the patch boundary-touch mask. Omitted/zero PCM/PLM
and PPM policies remain unchanged.

MPI two-level propagation, application input/default/restart policy and EB
remain outstanding. This increment tests hydro, not coupled chemistry or
transport. Extra ghost construction occurs only for a positive PCM/PLM
coefficient; it is not a performance qualification.

## Completed checks and failure

GNU Fortran 13.3.0, Ninja Debug, Python 3.12.3; reused
`build/spray-amr-reactive-debug`. Final focused CTest **1/1 passed**, 41.21 s
(41.26 s elapsed). Compiler warnings remain; no skips.

The existing AMR regression now checks coefficient 0.037 with PCM/PLM and
periodic/no-slip-plus-outflow boundaries, composite inventory versus boundary
flux after reflux (existing relative 2e-11 scale), a changed result relative
to zero, and negative-coefficient rollback of both levels and their
temperatures. The fine patch in this new fixture is interior. Existing entropy
wave, average-down, species closure, physical-boundary and API checks also run.
This is bounded conservation/API evidence, not external PeleC parity or
physical validation.

The first run failed at the pre-existing `public PLM ghost local-index
rejection` assertion after the new viscosity checks had passed. It supplied
index -2, which the existing four-layer PPM implementation already accepts
(valid minimum -3, confirmed in parent source). The test now supplies -4,
outside that existing range. No production validation or tolerance was relaxed.
The corrected regression passed on rerun.

Mechanism provenance, Python syntax, project contract and whitespace checks
passed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_amr_reactive_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^regression_amr_reactive_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/pcm-plm-two-level-viscosity-checks/` contains `build.log`,
`rebuild.log`, failed `LastTest-first.log` and final `LastTest.log`.
