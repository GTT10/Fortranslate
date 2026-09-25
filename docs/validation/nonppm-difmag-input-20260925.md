# PCM/PLM coefficient input and checkpoint identity — 2026-09-25

## Source and scope

Code `7a173ba`, parent `2834445`. Shared regular gas configuration accepts
explicit nonnegative coefficients for PCM/PLM in applications that enable
`allow_difmag`. Automatic -1 remains 0.1 for PPM/pelec and zero otherwise.
The existing application drivers already forward resolved coefficients.
Regular spray applications have separate namelists and remain outstanding.
Uniform reference-default policy and EB dispatch also remain outstanding.

Two-level checkpoints append the existing `difmag-v1` resolved coefficient
identity for positive PCM/PLM settings. Zero settings retain the legacy plain
reconstruction record, including automatic -1 and explicit zero. PPM identity
is unchanged. Tree checkpoints canonicalize zero PCM/PLM to the legacy -1
sentinel; positive values remain captured in the configuration context. Tree
context construction now rejects nonfinite coefficients and negative values
other than -1. Checkpoint headers and state payloads are unchanged.

## Completed checks

GNU Fortran 13.3.0, Ninja Debug, Python 3.12.3; reused fixed tests-enabled
`build/spray-amr-reactive-debug`. Focused CTest **4/4 passed**, total 41.09 s:
regular input 40.77 s, legacy two-level checkpoint 0.12 s, tree checkpoint
0.05 s, positive PLM two-level checkpoint 0.11 s. No failures or skips.
Compiler warnings remain.

The regular executable runs PCM/PLM with automatic, explicit zero, -1 and
positive 0.037 inputs: automatic/zero/-1 CSVs match exactly, and positive
output differs. Existing PPM input/default checks also pass. Two-level
checkpoint tests check positive PLM round-trip, rejection of coefficient 0.05
against 0.037 without changing state, and legacy automatic versus zero
identity. Tree context tests check zero equivalence and positive distinction.
The rest of those existing checkpoint tests run too.

These are fixed regular application and checkpoint-library checks. Selected
executables, changed-rank application restart and coupled particle runtime
were not independently rerun here. No external solver or physical validation
claim. Mechanism provenance, Python syntax, project contract and whitespace
checks passed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_3d test_amr_reactive_3d_checkpoint test_amr_patch_tree_checkpoint_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(reactive_3d_difmag_input|unit_amr_reactive_3d_(difmag_)?checkpoint|unit_amr_patch_tree_checkpoint_3d)$' --output-on-failure
python3 -m py_compile tools/check_reactive_difmag_input.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/nonppm-difmag-input-checks/build.log`, copied `LastTest.log`,
and input fixture paths recorded in that test log.
