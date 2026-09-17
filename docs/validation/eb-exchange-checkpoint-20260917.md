# EB checkpoint exchange vectors — 2026-09-17

## Change and boundary

Implementation: `7398e278630232eb2d688b9378b92fbcd46fc4e8`, based on
`4936b4e`. The checkpoint reader/writer accept optional complete boundary,
wall and chemistry exchange vectors. With all three supplied, fixed,
selected legacy-element and selected generic-element records use schemas
4, 5 and 6 respectively. Without these arguments, schemas 1/2/3 retain their
layout. Readers require the schema matching the requested context; no
implicit upgrade or migration is performed.

The `EXCHANGES` block precedes `STATE`. All vectors must have one entry per
conserved component and finite values. Writer preflight rejects incomplete,
wrong-sized or nonfinite data before opening the destination. Readers keep
all candidate data private until the entire record, state/EOS and end marker
have passed, then publish state and exchanges together. This is record
validation, not a proof that arbitrary finite ledgers balance the state.
Storage-failure atomic replacement is not added.

Single-level application accumulation, exchange-corrected component/element
gates and the earlier density-contrast methanol failure remain open. This
change does not modify solver updates or connect these optional arguments in
applications.

## Independently completed checks

GNU Fortran 13.3.0, `/usr/bin/gfortran`, Debug, Ninja, existing
`build/selected-tree-mpi-debug` configuration. No toolchain settings changed.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_reactive_eb_3d_checkpoint --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_reactive_eb_3d_checkpoint$' --output-on-failure --output-log build/eb-exchange-checkpoint-checks/ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

- Focused checkpoint unit: **1/1 passed**, 0.14 s (total 0.17 s),
  zero failures or skips. Includes existing legacy and borrowed-unit checks.
- Added schemas 4/5/6: exact signed nonzero vector and state round trips;
  invalid exchange marker, full-length nonfinite vector, truncated block and
  invalid temperature after exchange parsing reject without modifying targets.
- Incomplete and wrong-sized argument contexts reject. Nonfinite and invalid
  writer inputs preserve the existing checkpoint.
- Mechanism provenance, project contract and diff whitespace checks passed.
  No Python files changed; Python syntax checking is not applicable.
- Build output: `build/eb-exchange-checkpoint-checks/build.log`; test output:
  `build/eb-exchange-checkpoint-checks/ctest.log`. Build has real-comparison
  warnings; it completed successfully.

This is focused serialization/failure testing, not an application regression,
full configuration matrix, external PeleC parity or physical validation.
