# Single-domain EB viscosity input — 2026-09-25

## Source and scope

Code `eb5987b`, parent `d948f3c`. The fixed and selected single-domain EB gas
frontends accept finite nonnegative `difmag` in `&reactive_eb_3d`. The driver
forwards it through full/Strang composition into redistributed Euler hydro.
Default and explicit zero preserve the previous behavior. Negative values,
including -1, are invalid; this EB input has no automatic sentinel.

Positive coefficients extend the checkpoint solver identity with `difmag-v1`
and the coefficient. Zero retains the legacy solver identity and payload.
Changing a positive coefficient rejects restart. Applications sharing the
reader must explicitly enable this input; EB spray has not enabled it.
Tree/spray application wiring, uniform reference defaults and coarse/fine
geometry fidelity remain outstanding.

## Completed checks

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug; reused
`build/spray-amr-reactive-debug`. Focused CTests **2/2 passed**: existing
checkpoint unit 0.14 s and input/restart 4.60 s (4.78 s total). No failures or
skips. Compiler warnings remain. Selected frontend was edited but not built
or run in this increment.

The fixed density-sheet/state-redistribution fixture checks coefficient 0.037
continuous versus interrupted/resumed CSV equality, changed 0.05 rejection,
negative/NaN rejection, omitted versus explicit-zero equality and legacy-zero
restart. It confirms that the interrupted output differs from the completed
reference. It does not independently compare positive versus zero dynamics.
This is bounded application/restart evidence, not external solver parity or
physical validation. Provenance, Python syntax, contract and diff checks passed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_eb_3d test_reactive_eb_3d_checkpoint --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(reactive_eb_3d_difmag_input|unit_reactive_eb_3d_checkpoint)$' --output-on-failure
python3 -m py_compile tools/check_eb_difmag_input.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-difmag-input-checks/build.log` and `LastTest.log`.
Fixture outputs: `build/spray-amr-reactive-debug/tests/reactive_eb_3d_difmag_input/run-9eb6xygw/`.
