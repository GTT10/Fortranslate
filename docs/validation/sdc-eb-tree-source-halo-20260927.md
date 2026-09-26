# EB tree SDC source halos — 2026-09-27

## Source and behavior

Implementation `04c831a4`, based on `7787aac3` on local `main`.
PeleC pin remains `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

`fill_sdc_eb_tree_source_halo_3d` fills signed packed source fields using the
EB tree integer layout. Same-level peers take priority over PCM ancestor
donors. Optional ancestor start/end fields and fractions interpolate only
ancestor donors; peers and local cells retain current rates. The field has no
EOS or species-positivity constraint. This stage does not average down sources.

Periodic axes wrap; other accepted physical boundary labels use nearest-cell
scalar extrapolation, without momentum reflection or prescribed inlet state.
Widths 1 through 6 are accepted. Covered donor cells produce zero without
reading their rate. An optional target volume-fraction mask must match the
central patch and rejects an open target whose donor is covered. The complete
halo publishes only after all referenced values pass finite checks.

## Completed focused verification

GNU Fortran 13.3.0, Ninja Debug, fixed tests-enabled
`build/spray-amr-reactive-debug`. Existing compiler warnings remain.

- Mechanism provenance preflight passed for `mechanisms/h2o2_full.json`.
- `test_eb_tree_halo_3d` built successfully.
- CTest `unit_eb_tree_halo_3d`: 1/1 passed, 0.26 s test / 0.30 s total.
- New source tests: three levels/five patches, two signed components,
  three halo layers, exact cross-parent peer / parent / grandparent donors,
  time fractions 0, 0.5, 1, periodic corners, and scalar extrapolation for
  slip, no-slip, isothermal-wall and inflow labels.
- Solid-cell NaNs are ignored. Inconsistent open target/covered donor and a
  late nonfinite fluid donor both reject without altering the output halo.
- Existing gas halo checks remain in the same executable.
- Python syntax, project contract and `git diff --check` passed.

Commands from repository root:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_halo_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_halo_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Logs: `build/sdc-eb-tree-source-halo-checks/{provenance,build,ctest}.log`.
No selected build, MPI run, full matrix, installed-binary audit or physical
validation was performed. This is a source-halo interface, not a coupled tree
SDC advance. The existing gas halo still averages its snapshot by default;
the provisional SDC hydro stage must explicitly address that behavior when
connected. The scalar boundary and PCM interpolation choices match the current
EB/SDC interfaces; this record does not establish full pinned PeleC boundary
or coarse/fine interpolation parity.
