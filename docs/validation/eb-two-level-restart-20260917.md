# Two-level EB application checkpoint continuation — 2026-09-17

## Source and scope

Implementation: `37bb63646078b9ad4df82303367c9b420749eb6c`, based on
`6ba7a76` on local `main`. This connects checkpoint continuation to the
fixed-elementary, serial, one-fine-patch application. It does not qualify
selected mechanisms, MPI EB, multiple patches, subcycling or physical accuracy.

The V1 composite file records topology, chemical integrator, fixed maximum
step, time, step count, minimum accepted step/transport limiter, initial
composite inventories and cumulative boundary/wall/chemical changes. Two
embedded records retain existing mechanism, thermodynamic, transport, geometry
and EOS checks. Their diagnostic inventories describe each saved patch;
the outer record owns the composite initial state and cumulative budgets.

The reader stages both patches, requires the closing marker and EOF, checks
snapshot inventories, composite balance and parent averages, then publishes
state and diagnostics together. Final-time and maximum-step limits may extend;
physics, geometry, topology, integrator and fixed step must match. The writer
checks balance and serializes both patches to scratch before opening the
requested destination. Replacement is not atomic against storage failures.

The checkpoint serializer explicitly accepts redistribution `none` for plain
Euler patches, including a plane outside a patch or aligned to its faces.
Existing redistributed-mode restrictions, numerical tolerances and the
single-level application parser are unchanged. Input/output alias checks are
lexical, as in the existing application; they do not resolve filesystem links.

## Independently executed checks

Build: `build/selected-tree-mpi-debug`, Ninja, GNU Fortran 13.3.0, Debug;
configured MPI wrapper `/usr/bin/mpifort` (Open MPI 4.1.6). This application/test
is serial. Build and CTest ran with `LD_LIBRARY_PATH` and `LIBRARY_PATH` unset.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug \
  --target pelef_reactive_eb_two_level_3d test_reactive_eb_3d_checkpoint --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug \
  -R '^(regression_reactive_eb_two_level_3d|unit_reactive_eb_3d_checkpoint)$' --output-on-failure
python3 -m py_compile tools/check_reactive_eb_two_level_3d.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final focused run: **2/2 passed**, zero failures/skips, 1.98 s total
(checkpoint unit 0.07 s, application regression 1.88 s).
Provenance, Python syntax, project contract and whitespace checks passed.
Pinned mechanism source SHA256:
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

The application checker stops after one step at 1e-8 s, restores with extended
limits, and reaches three steps at 3e-8 s. Both final state CSV and complete
inventory CSV are byte-identical to uninterrupted execution. Normalized
inventory residual: `4.6057075715980637e-16`; minimum limiter: `1`.
It also rejects second-record truncation, altered cumulative mass exchange,
trailing data, invalid terminator, changed time step/transport setting and
checkpoint/input aliasing. Existing output files remain intact on rejection.

Logs: `build/eb-two-level-restart-checks/` (`initial-build.log`, `build.log`,
`final-build.log`, `ctest.log`, `final-ctest.log`, `verified-ctest.log`).
The last file owns the counts above. Earlier successful runs preceded added
negative cases and the writer balance preflight; they are not additional tests.
Artifacts and per-invocation logs:
`build/selected-tree-mpi-debug/tests/reactive_eb_two_level_3d/`.
No full configuration matrix, installed-binary audit, external solver parity
or physical validation was rerun.
