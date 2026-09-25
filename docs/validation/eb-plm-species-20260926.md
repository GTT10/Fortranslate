# Independent EB PLM species weights — 2026-09-26

Implementation `6657dc57`, parent `17f75071`, local main, version 0.246.0.

The EB normal predictor now disables the shared PLM helper's species
renormalization. Pinned PeleC `PLM.H::pc_plm_d_eb_passive` reconstructs each
species independently at open faces; its limited slopes need not sum to zero.
The downstream primitive Riemann adapter accepts these weights. Final CTU
species-flux normalization remains in place, including zero artificial viscosity.

The new regression uses three positive species with normalized cell values
`(0.1,0.3,0.6)`, `(0.2,0.2,0.6)`, `(0.4,0.1,0.5)`. The pinned order-two limiter
gives slopes `(1/5,-2/15,0)`, so zero-time lower/upper face sums are `29/30` and
`31/30`. Explicit expected component values check three axes and all four
open/closed-face combinations. Closed outputs remain zero.

## Verification

The added test failed before the repair at `lower EB species retain independent
slopes` (0.88 s test / 0.93 s total). After repair,
`regression_reactive_eb_hydro_3d` passed **1/1**, 23.37 s test / 23.41 s total,
with no skips. This includes existing CTU, redistribution, conservation and
rollback checks. Provenance, Python syntax, project contract and whitespace
checks passed. No tolerance or golden data changed.

GNU Fortran 13.3.0, Python 3.12.3, existing tests-enabled Ninja Debug build
`build/spray-amr-reactive-debug`; existing compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_eb_hydro_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^regression_reactive_eb_hydro_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-plm-species-checks/{provenance.log,build-before.log,ctest-before.log,build.log,ctest.log,regression.log}`.
Reference inspected using `git show` at pinned PeleC commit
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` in
`/mnt/sn850/Projects/gtt2/external/PeleC`.

This is a source-grounded focused regression, not a compiled external EB oracle
or full application/physical validation. MPI, selected mechanisms and installed
binaries were not rerun for this change.
