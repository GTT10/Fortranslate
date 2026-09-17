# Multiple-patch EB application — 2026-09-17

Implementation: `cfcb0eeff5b28277a24a39feeb327ed706b5b52a`, based on
`51cda2a` on local `main`.

`pelef_reactive_eb_two_level_3d` now accepts 2–16 static nonoverlapping fine
boxes through `patch_count`, `patch_lower` and `patch_upper`. It connects the
existing multiple-patch reaction/transport/hydro transaction, includes all
patches in automatic interval estimates and composite budgets, and adds a
patch-index column to multiple-patch state CSVs. The original one-patch path
and its input/output schema remain available.

Multiple-patch checkpoint V1/V2 bind fixed/automatic common stepping and store
all patch records, ordered topology and accumulated diagnostics together.
Restore validates every record, composite balance and parent averages before
publishing state. Patch order and boxes must match. Writes stage validation in
a scratch stream; replacement is not atomic against storage failure.

Scope is fixed elementary chemistry, serial, static two levels and a common
interval. Multiple-patch subcycling is explicitly rejected. Selected-mechanism
application support, dynamic/multilevel EB, MPI, general joined/split transverse
viscous-stencil parity and physical validation remain outside this result.

## Independently executed checks

GNU Fortran 13.3.0, Debug/Ninja, existing `build/selected-tree-mpi-debug`
(configured with `/usr/bin/mpifort`, OpenMPI 4.1.6). This regression is serial.

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target pelef_reactive_eb_two_level_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^regression_reactive_eb_two_level_3d$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/check_reactive_eb_two_level_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, zero failures/skips, 7.65 s (7.68 s total).
Provenance, Python syntax, contract and whitespace checks passed.
Mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

The added fixture splits the existing refined box into two adjacent patches
at ratio 2, with automatic interval selection, explicit reaction and transport.
It completes three steps to 3e-8 s without rejected trials, with normalized
composite inventory residual **2.3501078911122276e-16** and minimum transport
limiter 1. The 104 active-cell CSV rows have unique centroids, finite values,
positive density/temperature and species closure. Independent CSV integration
matches the reported final inventories.

Stopping after one step and resuming produces byte-identical state and budget
CSVs to the uninterrupted run. A truncated second-fine-patch record and changed
extra-patch topology reject restart while preserving existing outputs;
overlapping boxes reject initialization. Existing one-patch checkpoint,
automatic selection, retry and subcycling cases remain in the same regression.
The public example is `cases/reactive_eb_two_level_3d/multipatch.nml`; the
checker generates its own equivalent fixture with temporary output paths.

Logs: `build/eb-multipatch-app-checks/{build,final-build,ctest}.log`.
Generated inputs, restart files and CSVs:
`build/selected-tree-mpi-debug/tests/reactive_eb_two_level_3d/`.
No broader matrix, installed-binary audit, external solver parity, convergence
study or physical validation was rerun.
