# EB application reference viscosity default — 2026-09-25

Code `8ea199f`, parent `f944d28`. Shared single-domain EB configuration and
EB tree/two-level namelists default `difmag` to 0.1. This covers gas and spray
application inputs, including shared selected/MPI tree source. Values must be
finite and nonnegative; -1 remains invalid in EB inputs. Explicit zero retains
the earlier zero correction and checkpoint identity. Automatic/default restart
from a zero-coefficient checkpoint now rejects instead of changing dynamics.

The shared reader accepts positive values by default now that its application
callers support the coefficient. An explicit `allow_difmag=.false.` still
rejects them. Direct-library omitted arguments remain on the earlier zero policy
pending unification; this change is at the application/configuration boundary.

## Completed checks

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug; reused
`build/spray-amr-reactive-debug`. Fixed single-domain/tree/two-level gas and
serial/MPI EB spray executables built. **6/6 CTests passed** in two invocations:
single gas 6.49 s, tree 37.27 s, two-level 5.66 s, multipatch 6.17 s, serial
spray force/restart 1.10 s (56.73 s total); configuration unit 0.00 s (0.05 s
elapsed). No failures or skips. Compiler warnings remain. Provenance, Python
syntax, contract and whitespace checks passed.

Gas input checks confirm omitted versus explicit 0.1 exact state CSV equality,
restart across those equivalent settings, explicit-zero continuation, and
zero-to-default rejection. Existing 0.037 continuation, changed-coefficient,
negative and NaN checks pass. Serial spray runs the new default with exact
gas/particle/film restart checks. The last reader opt-in-default adjustment
was checked by rebuilding/running the configuration unit; application callers
already explicitly enable the coefficient and were not rerun after that edit.

MPI runtime, selected builds, historical binary checkpoints, geometry fidelity
and the full regression matrix were not checked. No full external solver parity
or physical validation is claimed. Historical zero-policy fixtures must request
zero explicitly; no golden hashes or tolerances were changed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_eb_3d pelef_reactive_eb_tree_3d pelef_reactive_eb_two_level_3d pelef_eb_spray_3d pelef_mpi_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(reactive_eb_3d_difmag_input|reactive_eb_tree_3d_difmag_input|reactive_eb_two_level_difmag_(reference|multipatch)|spray_eb_force_application)$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_reactive_eb_3d_config --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^unit_reactive_eb_3d_config$' --output-on-failure
python3 -m py_compile tools/check_eb_difmag_input.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-default-difmag-checks/` contains `build.log`,
`config-build.log`, `LastTest-app.log` and `LastTest-config.log`.
