# MPI recursive dynamic-tree integration — 2026-09-17

## Source and purpose

Base `45974f5be639c6334cabebd69c2a8ea8aeee8507`, local `main`.
Extends the existing scheduled-restart and automatic-adaptation checkers with
recursive mode. This is an integration qualification change; solver algorithms,
mechanism inputs, tolerances and golden hashes are unchanged. Dynamic tests
use hydro only to isolate geometry, interval scheduling and restart; the static
full reaction/transport/hydro application was checked separately.

SHA-256:

- `tools/check_mpi_amr_tree_restart_3d.py`:
  `12f10bf35d62d618d7490355c79d1e4d5185f50e4bd8ad4fd3617eb359f6a267`
- `tools/check_mpi_amr_tree_adaptation_3d.py`:
  `504f5042c52d1cb55e941086c21fc1b64db97c2509d6584b080d5325647f05fd`

## Commands

GNU Fortran 13.3.0 Debug / OpenMPI 4.1.6, existing
`build/selected-tree-mpi-debug`; no compiler/backend changes.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target pelef_mpi_amr_tree_reactive_3d pelef_amr_tree_reactive_3d test_selected_elementary_mpi_amr_tree_3d test_selected_elementary_amr_tree_3d -j 4
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_mpi_amr_tree_subcycle_(schedule|adaptation)_3d$' --output-on-failure
python3 -m py_compile tools/check_mpi_amr_tree_restart_3d.py tools/check_mpi_amr_tree_adaptation_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

## Completed results

Build targets passed. CTest **2/2 passed**, 75.43 s total: scheduled 33.25 s,
automatic 42.15 s; no failures or skips. Provenance, Python syntax, project
contract and whitespace checks passed. Both automatic runs completed three
coarse steps with ten nodes at three levels; maximum reported remap error was
1.0569e-16. Scheduled maximum remap error was 1.1179e-16.

## Coverage

Scheduled PCM execution uses the existing three-event case: ratio-changing
remap at 2.5e-10 s, root-only coarsening at 5e-10 s and branched refinement at
1e-9 s. Three coarse intervals end exactly at those events. Fixed and selected
frontends must match serial CSV/budget bytes and step counts, including 2-to-3
rank restart, serial interchange and exact checkpoint bytes. A checkpoint at
the root-only event exercises empty ranks on reload. Mode/schedule/physics,
truncation/ledger and path-alias guards remain checked. The final refinement
occurs at final time; there is no subsequent interval on that final tree.

Automatic temperature adaptation uses characteristic PLM/MC, an off-center
hotspot and two refinement levels with ratios 2,2. It starts from the root,
checks tagging every coarse step and verifies final exposed-cell coverage,
level two, exact serial CSV/budget bytes and step count, 2-to-3 restart,
serial interchange and checkpoint-byte identity. Adaptation-counter tampering,
changed mode/threshold/cadence and node-budget exhaustion reject while preserving
existing outputs. Fixed and selected final output bytes must also match.

Logs: `build/amr-tree-mpi-subcycle-dynamic-checks/{build,ctest}.log`.
Artifacts: `build/selected-tree-mpi-debug/tests/mpi_amr_tree_subcycle_{schedule,adaptation}_3d/`.

## Limits

This does not qualify dynamic recursive molecular transport/chemistry, inlet
bound selection, general temporal convergence or production scaling. No full
matrix, installed audit, external-solver parity or physical validation was run.
