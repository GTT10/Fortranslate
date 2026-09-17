# MPI AMR reacting spray boundary composition — 2026-09-18

Implementation: `84a11a94e64a81893909f094aa254cd533ee1933`, based on `c52cd88`.

The reacting wrapper now accepts optional face types, a conserved boundary
ledger and existing gas inflow/thermal/moving-wall data. Default faces remain
periodic. Both particle halves and the gas split receive the same faces; their
integrated inventory changes sum into a private ledger. State, diagnostics,
IDs and the optional ledger publish together after all stages succeed. Failed
calls preserve the incoming ledger. Shape and overlong face checks are
collective; the particle interval checks face validity and rank agreement.

This connects the library composition. The application and checkpoint still
use periodic boundaries. Wall faces remain specular for particles; they do not
introduce films, wall heat exchange with droplets or EB impacts.

## Verification

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja, `/usr/bin/gfortran`,
`/usr/bin/mpifort`, `/usr/bin/mpiexec`. Reused unchanged fixed tests-enabled and
selected H2/O2 tests-disabled configurations.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_mpi_amr_reactive --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_mpi_amr_reactive_(2|4)$' --output-on-failure --output-log build/spray-amr-reactive-boundary-checks/ctest-final.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final focused CTest: **2/2 passed**, no failures or skips, two ranks **57.74 s**,
four ranks **55.98 s**, total **113.76 s**. Fixed and selected library builds,
provenance preflight, Python syntax, project contract and whitespace checks
passed. Existing compiler warnings remain.

The three-level fixture retains periodic common/recursive conservation,
implicit chemistry continuation and gas-stage failure rollback, now including
a sentinel ledger. Added slip-wall/outflow common/recursive cases retain
active drag, heat, evaporation and gas viscosity/conduction/diffusion. Composite
gas-plus-liquid inventory differences match the combined boundary ledger at
3e-12 of inventory/ledger scale, with an exiting parcel removed. The test timeout
increased from 60 to 120 seconds to accommodate the added coupled cases.

The first two tests failed because initially zero transverse momentum was
compared against a roundoff-sized final value (about 2e-23 kg m/s); a diagnostic
rerun found absolute residuals below 1.1e-25 kg m/s in those components. The
fixture now retains nonzero transverse particle momentum, so each relative
budget assertion has a physical scale. No tolerance or numerical implementation
was changed to address that failure. Logs retain the initial failed run and
diagnostic output separately from final results.

Selected runtime compilation is compile evidence only. Prescribed inflow and
thermal/moving-wall arguments are forwarded but not newly exercised here;
application/restart integration, external parity and physical validation remain
open. No full matrix was run.

Artifacts: `build/spray-amr-reactive-boundary-checks/` contains `build.log`,
`build-selected.log`, `build-diagnostic.log`, `ctest.log` (initial failures),
`diagnostic.log` and `ctest-final.log`. Fixed executable:
`build/spray-amr-reactive-debug/tests/test_spray_mpi_amr_reactive`.
