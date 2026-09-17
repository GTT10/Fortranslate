# Single-level EB application exchange budgets — 2026-09-17

Implementation: `20ca792a666f49a8c033551ef186002a4fe7b830`, based on
`698a0b4`. No numerical update kernels or mechanism data changed.

## Behavior

The shared single-level application accumulates accepted hydro boundary,
wall and chemical exchange vectors. Component residuals are final minus
initial minus all three exchanges. Element residuals subtract external
boundary/wall exchange only, retaining the chemical elemental-conservation
check. A signed amount-to-element projection supports generated composition
metadata and the fixed H/O/N mechanism.

All components are now checked, including wall-normal momenta and reacting
species formerly skipped. Relative gates remain 5e-11 for components and
5e-10 for elements. Component normalization includes initial/final L1 amounts
and magnitudes of exchanges. The newly checked species budgets have an
absolute summation allowance of 64 epsilon times accepted steps times total
initial absolute species mass; the existing generic-element roundoff allowance
is retained. Raw residuals are printed separately from these corrected errors.

Fixed and selected application checkpoints now use schemas 4 and 6. Older
application records are rejected; library callers may still use legacy schemas
without exchange arguments. On restart, the restored budget is checked before
stepping or writing output; checkpoint writes also require the current budget
to pass. The library reader remains transactional for malformed records;
application budget rejection aborts before external output. This does not add
atomic file replacement on storage failure or cryptographic ledger integrity.

The fixed restart checker retains the existing schema-1 payload golden SHA:
it checks the new exchange structure and strips only that added block plus the
schema digit before comparing the original bytes. The installed full-selected
checker similarly accounts for added source element and exchange metadata,
with explicit element-table checks. No golden hashes were replaced.

## Density-contrast result

The methanol fixture now uses regular/cut density 0.31/0.93 kg/m3, the same
contrast that previously failed the zero-exchange mass gate. Five steps reach
1e-11 s with the implicit 21-species, 124-reaction FFCM1_Red bundle.

| Quantity | Result |
| --- | --- |
| Initial mass | 2.4800000000000000e-13 kg |
| Final mass | 2.4799999998720304e-13 kg |
| Integrated boundary mass exchange | -1.2797128019602285e-23 kg |
| Wall / chemical mass exchange | 0 / 0 kg |
| Raw mass budget residual | 1.6330623461429972e-28 kg |
| Maximum component budget error | 6.5849288150927309e-16 |
| Element error after existing roundoff allowance | 0 |

The measured boundary exchange explains the previous mass gate failure.
Uninterrupted and checkpoint/resumed runs have exactly equal final CSV bytes
and printed cumulative budgets. Corrupting checkpoint boundary mass to 1e-4
is rejected before continuation, preserving the prior output and checkpoint.

## Completed verification

GNU Fortran 13.3.0 (`/usr/bin/gfortran`), Debug/Ninja, existing
`build/selected-tree-mpi-debug`. Optional methanol bundle:
`build/port-methanol-reference/ffcm1_red.json` (source SHA
`ddbba8c97d929530cfb213da45485cca4f979dccf97c1063721bd855d7b3f259`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_selected_methanol_eb_3d pelef_reactive_eb_3d test_selected_elementary_reactive_eb_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^(regression_selected_methanol_eb_3d|regression_reactive_eb_restart_3d_(reference|checkpoint_stop|restart|check)|regression_selected_reactive_eb_3d_restart_(reference|checkpoint_stop|restart|composition_mismatch_rejected|check))$' --output-on-failure --output-log build/eb-app-exchange-checks/ctest.log
python3 -m py_compile tools/check_selected_carbon_eb_3d.py tools/check_reactive_eb_3d_restart.py tools/check_selected_reactive_eb_3d_restart.py tools/check_installed_full_selected_reactive_eb_3d_restart.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final focused run: **10/10 passed**, 14.40 s, zero failures/skips. Includes
fixed and selected H2/O2 continuation and selected composition rejection,
plus methanol exact state/budget continuation, element-metadata rejection and
finite-but-inconsistent budget rejection. Provenance, Python syntax, contract
and diff checks passed. An initial build failed on one overlong Fortran line;
that line was wrapped and the final build passed.

Build/test logs: `build/eb-app-exchange-checks/`. Carbon artifacts:
`build/selected-tree-mpi-debug/tests/selected_methanol_eb_3d/`.
The installed full-selected checker was syntax-checked, not exercised against
a fresh install. Release, the full configuration matrix, general geometry
regressions, external PeleC parity and physical/accuracy validation were not
rerun. General EB multilevel/dynamic/MPI integration remains broader open work.
