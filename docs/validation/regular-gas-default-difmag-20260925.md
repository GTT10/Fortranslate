# Regular gas reference viscosity default — 2026-09-25

Code `5812771`, parent `a2495cc`. Pinned PeleC
`bf0e1fd15040f0f5609cd9042b9f1b868e0e95f8` defines `difmag=0.1` in
`Source/Params/_cpp_parameters` and generated `pelec_defaults.H`, independent
of reconstruction. The shared regular gas resolver now maps omitted/-1 to
0.1 for every solver/reconstruction. Explicit nonnegative values override it.
The single-domain application now uses this resolver rather than duplicating
a PPM/pelec-only policy. Serial/MPI AMR applications already call the resolver.

This deliberately changes previously zero automatic PCM/PLM and non-pelec
application behavior. Set `difmag=0` to request the former zero correction.
Checkpoint identity uses resolved coefficients: automatic and explicit 0.1
are equivalent, while zero differs. Older zero-policy checkpoints require
explicit zero, rather than silently continuing under the new default.
No golden hashes or tolerances were changed.

## Completed checks

GNU Fortran 13.3.0, Python 3.12.3, Ninja Debug; reused
`build/spray-amr-reactive-debug`. **4/4 CTests passed**, 46.83 s total:
regular input 46.53 s, two-level checkpoint 0.11 s, tree context 0.05 s,
explicit-coefficient checkpoint 0.11 s. No failures or skips. Provenance,
Python syntax, project contract and diff checks passed. Compiler warnings remain.

Fixed regular PeleC-solver input checks PCM/PLM/PPM omitted, explicit 0.1 and
-1 exact CSV equivalence; explicit zero differs. Explicit 0.037 remains valid.
Two-level checkpoint tests reject changed coefficients with state preserved.
Tree context tests distinguish zero from automatic and equate explicit 0.1.
Selected and MPI application runtimes, historical checkpoint files and the
full regression matrix were not independently run. Existing historical
fixtures that depend on zero must request zero explicitly; no claim that the
entire matrix passes is made.

Direct library omitted coefficients, separate regular spray defaults and EB
defaults retain their previous policies pending further unification. MPI AMR
spray already uses the shared resolver. Reference geometry fidelity and full
external solver parity remain open.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target pelef_reactive_3d test_amr_patch_tree_checkpoint_3d test_amr_reactive_3d_checkpoint --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^(reactive_3d_difmag_input|unit_amr_patch_tree_checkpoint_3d|unit_amr_reactive_3d_checkpoint|unit_amr_reactive_3d_difmag_checkpoint)$' --output-on-failure
python3 -m py_compile tools/check_reactive_difmag_input.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/regular-gas-default-difmag-checks/build.log`, `LastTest.log`.
