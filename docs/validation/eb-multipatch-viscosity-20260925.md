# Multipatch EB viscosity — 2026-09-25

Code `ded76ea`, parent `840b87b`. Multipatch Euler and full reaction/transport
composition accept optional finite nonnegative `difmag` (default zero).
Validation precedes private chemistry work. Hydro applies the coefficient on
coarse and fine updates before flux accumulation, peer reconciliation and reflux.
Common and recursive steps request complete patch halos, including temporally
interpolated parent donors, then overlay same-level peer states/volume fractions.
Molecular transport does not use this coefficient.

The optional volume-fraction argument of peer overlay covers every ghost slot
inside the peer box, including edges/corners and covered cells. Covered peer
states become zero without EOS conversion. Existing shared-face metric checks
remain. Halo outputs commit together. Omitting the argument retains face-only
behavior. Parent geometry is still PCM; fine reconstruction and mixed physical
boundary/corner fidelity against PeleC are not established.

## Completed checks

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug; reused
`build/spray-amr-reactive-debug`. **3/3 CTests passed** across two invocations:
existing hydro/transport fixture 11.99 s, positive coefficient fixture 12.33 s
(24.37 s elapsed), positive full composition 12.83 s (12.87 s elapsed).
No failures or skips. Compiler warnings remain. Provenance, Python syntax,
project contract and whitespace checks passed.

Coefficient 0.037 covers common/recursive stepping, ratios 2/3, full-face,
partial-face and corner-touching layouts. Existing inventory bounds and
split/joined patch comparisons (3e-12 scaled bounds) pass. The full fixture
covers explicit/implicit chemistry, transport, exact patch-order invariance,
disabled-operator hydro parity and rollback after chemistry/transport halves.
Existing element conservation bounds pass without tolerance changes.

This is bounded library evidence. Older application coefficient inputs,
uniform defaults and reference geometry fidelity remain open. No full matrix,
external solver parity or physical validation is claimed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_multipatch_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_multipatch(_viscosity)?_3d$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_multipatch_full_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_multipatch_viscosity_full_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-multipatch-viscosity-checks/` contains `build.log`,
`full-build.log`, `LastTest-hydro.log` and `LastTest-full.log`.
