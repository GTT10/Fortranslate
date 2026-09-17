# Serial EB mixed-particle coupling — 2026-09-18

Implementation: `a1d923842ad4ac349db1cac0ecb06f3963c5adde`, based on `7f4f415`.

The interval connects EB sampling, multicomponent drag/heating/evaporation,
representative-drop ETAB/KH-RT breakup, specular drift and conservative gas
deposition. Common substeps use drag, component-loss, temperature-change and
cell-transit bounds. Rejected trials halve up to 30 times. All state remains
private until the complete interval succeeds. Wall momentum/energy enters a
separate particle-sign ledger. Failed intervals preserve state and ledger and
return zero accepted substeps.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; existing fixed tests-enabled and
selected H2/O2 tests-disabled configurations reused. Compiler
`/usr/bin/gfortran`, MPI wrapper `/usr/bin/mpifort`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_step --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_(step|field|drift|intersection)$' --output-on-failure --output-log build/spray-eb-step-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_mechanism_provenance.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **4/4 passed**, total **0.08 s**: step **0.02 s**, field **0.01 s**,
drift/intersection each **0.00 s**. Selected runtime compilation, provenance,
Python syntax, contract and whitespace checks passed. The initial test build
failed because the transport-loader call supplied an extra species argument;
correcting the call fixed compilation. No runtime failures or skips; no
tolerance changes.

The new fixture uses a 4 mm cube of 1 mm cells, a planar cut at x=1.25 mm,
1000 K N2 gas with nonzero momentum, and a 350 K water/peroxide parcel with
100 represented droplets, 1e-12 kg each, initially 1 micrometre from the wall.
The interval is 20 microseconds. Assertions require multiple substeps, heating,
evaporation, drag, reflection, positive wall impulse, a fluid endpoint and
unchanged covered cells. Every integrated gas-plus-liquid component closes
against the wall ledger at relative 3e-12 with a 1e-30 scale floor. Limiting the
interval to one substep fails and preserves gas, temperature, composition,
position and ledger. This fixture disables breakup; its newly connected
representative-drop branch is not separately requalified here.

Artifacts: `build/spray-eb-step-checks/build.log` (initial compiler failure),
`build-final.log`, `build-selected.log`, `ctest.log`.

The method applies sources before drift at updated velocity, replacing the
free-space source kernel's trapezoidal position update. It is first order;
convergence is unqualified. Feedback uses pre-substep positions. KH children
are rejected, zero-mass records retained, and trial endpoints must remain in
the root domain. This serial interval does not advance gas chemistry/transport
or hydrodynamics, form wall films, compose exterior events, or migrate particles
across MPI/AMR ownership. No external parity or physical validation was run.
