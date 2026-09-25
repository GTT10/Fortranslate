# Serial recursive/coupled EB tree viscosity — 2026-09-25

## Source and scope

Code `16dffd2`, parent `2e06648`. Serial recursive EB hydro accepts an optional
finite nonnegative coefficient. Each level step requests complete state and
volume-fraction halos using the existing ancestor time fractions and passes
them into EB Euler flux correction before flux accumulation/reflux. Full
reaction/transport composition forwards the coefficient to common or recursive
hydro; transport operators retain their separate molecular-viscosity controls.
The full wrapper validates the coefficient before private reaction work.

Omitted/zero coefficients retain previous behavior. MPI and application
input/default/restart wiring remain outstanding. The existing PCM ancestor
volume-fraction approximation and first-active-axis inlet corner policy remain;
their fidelity to reference geometry/corners is not established by these tests.

## Completed checks and failure

GNU Fortran 13.3.0, Ninja Debug, Python 3.12.3; reused
`build/spray-amr-reactive-debug`. Focused CTest **1/1 passed**, 48.91 s
(48.95 s elapsed). Initial build failed on two overlong test call lines;
wrapping them fixed compilation. No runtime failures or skips. Compiler
warnings remain.

The existing five-node z-plane fixture uses coefficient 0.037 with Mach inflow,
recursive stepping, four chemistry/transport combinations and explicit/implicit
reaction integrators. Combined state, temperature and diagnostics match the
independently called split stages exactly. Composite inventory balance retains
the existing scaled 3e-14 bound. Failure after private reaction/transport and
invalid inlet controls retain state and diagnostics. No tolerances were changed.
Common-step full composition was compiled but not independently run in this
increment. This is bounded split/hierarchy evidence, not external PeleC parity,
a full matrix or physical validation.

Mechanism provenance, Python syntax, project contract and whitespace checks
passed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_subcycle_full_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_viscosity_subcycle_full_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-tree-viscosity-split-checks/` contains initial `build.log`,
successful `rebuild.log` and copied `LastTest.log`.
