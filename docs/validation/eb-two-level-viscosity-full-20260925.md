# Coupled two-level EB viscosity — 2026-09-25

Code `517a0fd`, parent `a20b487`. Two-level Strang/full composition accepts
optional finite nonnegative `difmag`. The shared split wrapper validates it
before private chemistry/transport, then forwards it to common or subcycled
hydro, including the chemistry/transport-disabled shortcut. Molecular
transport coefficients are unchanged. Omitted values retain zero behavior.

## Completed checks

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug; reused
`build/spray-amr-reactive-debug`. Focused CTest **1/1 passed**, 17.83 s
(17.87 s total), no failures or skips. Compiler warnings remain. Mechanism
provenance, Python syntax, project contract and diff checks passed.

Coefficient 0.037 covers common/recursive stepping, explicit/implicit chemistry,
full reaction/transport, disabled chemistry, disabled transport and pure hydro.
Disabled branches match their hydro/Strang references exactly. Existing
composite inventory bound 3e-12 and element bound 2e-10 pass. Covered storage
remains unchanged. Negative coefficient rejects without changing states,
temperatures or ledgers. Existing bad-integrator and post-reaction hydro
failure rollback checks run with the positive coefficient.

This is bounded coupled-library evidence. Multipatch peer overlay, older
application input wiring, reference geometry fidelity and uniform default
policy remain open. No external solver parity or physical validation is claimed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_two_level_full_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_two_level_viscosity_full_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-two-level-viscosity-full-checks/` contains `build.log`,
`rebuild.log` and `LastTest.log`.
