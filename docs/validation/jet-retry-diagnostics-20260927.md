# Jet continuation and rejection diagnostics — 2026-09-27

Implementation: rejection logging in `app/pelef_reactive_eb_tree_3d.F90` now
prints step, accepted time, attempted interval and retry number, plus processed
injection events for the spray application. MPI prints these fields on rank
zero only. No acceptance rule, retry policy or numerical update changed.

The preceding full jet failed with retries explicitly disabled. Its last
checkpoint is step 2000 at `0.000359 s`; the previous log cannot locate the
subsequent rejected step. This does not prove timestep size caused the failure.
A separate continuation in `build/decane-jet-retry-run` uses the exact original
binary SHA-256
`31b374ef2bdc0385bb464e0313a50baaa1548cd0a6692d01c7d3abc0c63bcdda`.
Only `maximum_step_retries` changes from 0 to 12 and checkpoint interval from
1000 to 100; the restart points to the preserved original step-2000 checkpoint.
Existing transactional half-interval retries are thereby enabled. Boundary
conditions, event/rate files and physics are unchanged. The original failure
artifacts remain intact. Checkpoint hashes and control differences are stored
in `restart-provenance.json`; binary identity and launch are in `execution.json`.

The continuation started at `2026-09-27T09:25:07.422768+00:00`, driver PID
604587, launcher 604588, solver 604591. The live solver was independently
observed at approximately one minute elapsed, 99.6% CPU. This is a running
continuation, not a completed or repaired full-jet result. It deliberately uses
the earlier binary, so it does not contain the new rejection logging.

## Focused diagnostic verification

GNU Fortran 13.3.0 / OpenMPI 4.1.6, existing heptane Debug MPI/spray build:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray pelef_reactive_eb_tree_3d_selected --parallel 2 > build/tree-rejection-diagnostic-build.log 2>&1
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

All pass. A copy of the checked-in heated-wall case with deliberately excessive
`external_forcing=1e12,0,0` was run on two ranks from
`build/tree-rejection-diagnostic` with the selected MPI spray executable.
It exits 1 as expected, publishes no final CSV, and prints exactly one context
record: step 1, time 0, interval `9.9999999999999995e-8 s`, retry 0, one
processed injection event. Its terminal error remains `Tree step retry limit
exhausted`. `fortran.log` preserves the output. This tests rejection reporting,
not physical validity of that forcing. No full matrix or Python source change.
