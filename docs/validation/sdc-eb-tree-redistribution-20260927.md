# Composite EB tree flux redistribution — 2026-09-27

Implementation `d9641a76`, based on `a95acd59` on local `main`.
PeleC reference pin: `bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8`.

The serial library now provides `redistribute_sdc_eb_tree_flux_rate_3d`.
This is a composite source-rate operator; it is not yet called by the SDC driver.
It extends the existing single-patch flux redistribution to face-connected
leaf cells across peer and coarse/fine boundaries. Hidden parent cells neither
donate nor receive; their input rates are retained. Solid output rates are zero.

For each cut leaf, the operator forms an average of its original rate and
neighbor rates using physical fluid volumes. With volume fraction kappa, its
new rate is kappa times its original rate plus (1-kappa) times the average.
The volume-integrated difference is distributed to fluid neighbors in proportion
to their volumes. All neighborhoods read the original rates; incoming excess
is accumulated separately. The same weights apply to every signed component.
On a uniform grid this reduces to the existing volume-fraction-weighted operator.

Neighbors are located using integer tree coordinates. Both faces must be open;
solid cells and cells replaced by children are excluded. A cut leaf with no
fluid neighbor is rejected. Outputs publish only after the full candidate is
finite. Malformed inputs and nonfinite fluid rates preserve the caller's output;
unused solid rate NaNs are ignored. Physical domain faces are closed for this
operator; periodic wrapping is not yet implemented.

## Completed focused checks

GNU Fortran 13.3.0, Ninja Debug, fixed tests-enabled build
`build/spray-amr-reactive-debug`; existing compiler warnings remain.

- Mechanism provenance and the new test target built successfully.
- CTest `unit_eb_tree_redistribution_3d`: 1/1 passed, 0.01 s test / 0.06 s total.
- Root-only redistribution agrees with the single-patch implementation.
  Two fine peer patches tiling the root agree with one equivalent fine grid;
  maximum absolute component comparison threshold is 2e-14.
- A partially refined root conserves both signed rate components. Separate
  coarse and fine impulses demonstrate transfers in both directions across
  the interface, with integrated discrepancy below 3e-14.
- Relative composite residuals for root-only, partially refined and fully
  refined layouts: 2.7755575615628914e-17, 1.0856786298201405e-15 and
  7.1969637401713039e-16 respectively.
- Nonfinite fluid input and absent fields reject without publishing output.
- Python syntax, project contract and `git diff --check` passed.

The first test invocation failed because its child patch descriptor incorrectly
specified the parent's dimensions. The fixture was corrected before the final
checks; no numerical tolerance was changed. No skipped or interrupted checks
were counted.

## Reproduction

Artifacts: `build/sdc-eb-tree-redist-checks/{provenance,build,ctest,detail,contract}.log`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug \
  --target test_eb_tree_redistribution_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug \
  -R '^unit_eb_tree_redistribution_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Remaining work: connect the operator to SDC stage assembly with consistent flux
synchronization and ledger accounting; add periodic connectivity, MPI owner-local
transfers and recursive-interval redistribution accounting; extend state
redistribution. This focused operator check establishes neither time-stepping
stability nor physical validation or upstream AMR redistribution parity.
No selected-mechanism, installed-binary, configuration matrix or scaling audit
was run in this increment.
