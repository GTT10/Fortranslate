# Periodic serial AMR Dynamic LES — 2026-09-25

Implementation `c5088497fee03a54149435a8a8c766a701376ff0`, based on `9feb99b`.

The existing serial AMR LES SSPRK2 stepper and timestep selector now dispatch
`model='dynamic'` to the supplied-halo dynamic operator. Each patch receives
filter-radius-plus-five donor layers. Existing peer/coarse-fine synchronization,
reflux, averaging and EOS recovery retain their transaction boundary. Dynamic
filter controls default to box/2/1/arithmetic; physical boundaries reject.
Static coefficients do not control the dynamic closure. The ghosted operator's
new optional diffusivity output is a conservative maximum over the padded
workspace, so it can yield smaller intervals than an interior-only estimate.

This is the periodic serial library path. MPI Dynamic AMR, coupled/application
dispatch, nonperiodic dynamic boundaries and LES subcycling remain unconnected.

## Completed verification

GNU Fortran 13.3.0, Ninja Debug, existing fixed tests-enabled and selected H2O2
tests-disabled builds. Fixed `test_les_amr_tree` and selected
`pelef_selected_spray_tree_runtime` built. Selected execution was not run.
One direct `test_les_amr_tree` invocation passed, including existing static
Smagorinsky/WALE/Vreman cases and new periodic Dynamic cases on a root-only and
a nested multipatch tree.

Root dynamic state/temperature agree with regular-grid SSPRK2 at existing
`2e-13` scaled state / `1e-9 K` tolerances; boundary ledger agrees within `1e-14`.
Composite conservation holds within `2e-12` scaled tolerance. The timestep bound
is conservative against the regular operator and scales linearly with CFL.
Excessive timesteps preserve all node states, temperatures and diagnostics.
Invalid CFL preserves selector outputs. No tolerances changed.

The first build failed on one overlong Fortran argument-list line; wrapping it
fixed the error, and the final build and test passed. Provenance, Python syntax,
project contract and whitespace checks passed. Earlier MPI launcher handles
2018710/2018713 were still present on this turn's initial check; no new MPI test
was started or reported successful. No full matrix, external solver comparison,
physical validation or scaling check was performed.

## Commands and artifacts

Logs: `build/les-dynamic-amr-checks/` (`provenance.log`, `build.log` initial failure,
`build-final.log`, `build-selected.log`, `test.log`, `contract.log`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_les_amr_tree --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH build/spray-amr-reactive-debug/tests/test_les_amr_tree
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```
