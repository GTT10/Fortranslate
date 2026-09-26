# Provisional EB tree SDC hydro — 2026-09-27

## Source and behavior

Implementation `1ee06fef`, based on `647f4ca6` on local `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`compute_sdc_eb_tree_hydro_3d` connects conserved predictor sources to the
existing EB CTU hydro stage on every patch. It prepares six gas halo layers,
three source layers, flattening and physical boundary indices using the tree's
integer layout and extended geometry. Outputs are packed conservative rates
per fluid volume and time-integrated face flux densities (without aperture).
All outputs publish together after every patch succeeds.

The gas halo now has an optional `restrict_snapshot` control, defaulting true
for existing callers. The provisional SDC stage passes false to retain parent
states during iteration. It performs no average-down, redistribution or reflux;
those operations belong to the coupled driver. The interface supports common
steps; no selected-level or ancestor-time interface is connected here.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance preflight passed for `mechanisms/h2o2_full.json`.
- Hydro and halo unit targets built successfully.
- CTest hydro/halo: 2/2 passed, 42.94 s total (hydro 42.63 s, halo 0.26 s).
- New hydro check: z-plane cut, three levels, five patches, zero predictor.
  Reconstructing state as old state plus dt times rate and then applying the
  existing flux synchronization reproduces the common-step hydro conserved
  state within 2e-12 times max(1, abs(reference)). This tests assembly and
  synchronization compatibility, not physical convergence.
- Nonzero force/work predictor changes deep-patch hydro rates. Solid-cell
  source NaNs are ignored. An active source NaN rejects with every prior rate
  and x/y/z flux array unchanged.
- After adding a direct provisional-parent assertion, rebuilt and reran the
  halo executable separately: passed. With restriction disabled, root boundary
  gas equals the original root donor exactly despite different refined states.
- Python syntax, project contract and `git diff --check` passed.

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_hydro_3d test_eb_tree_halo_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_(hydro|halo)_3d$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_halo_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH build/spray-amr-reactive-debug/tests/test_eb_tree_halo_3d
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Logs: `build/sdc-eb-tree-hydro-checks/` (`provenance.log`, `build.log`,
`ctest.log`, `halo-build.log`, `provisional-halo.log`).

No coupled EB tree SDC trajectory, recursive stepping, selected build, MPI run,
full configuration matrix, installed-binary audit or physical validation was
performed. Existing CTU and scalar source boundary policies are reused; these
checks do not establish complete pinned PeleC algorithm parity.
