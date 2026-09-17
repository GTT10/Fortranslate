# Multilevel EB checkpoint records — 2026-09-18

Implementation: `5b66385dc3fbaed9ccd6511ca7cd14c44e9dbbbc`, based on
`7870b41`. Core module: `eb_reactive_tree_checkpoint_3d_mod`.

## Contract and integration boundary

`reactive_eb_tree_checkpoint_record_3d` reads/writes one formatted stream
envelope. `eb_tree_checkpoint_3d` holds time, steps, minimum interval/limiter,
initial composite inventory/scale and cumulative boundary, wall and chemical
exchange. V1 stores exact node count and topology descriptors followed by
metadata and existing single-grid EB records for every node. Embedded records
reuse physics, geometry, EOS and optional selected-mechanism identity checks;
their inventory fields describe the stored node snapshot.

The caller supplies the expected initialized tree topology and geometry.
Reads stage the complete tree and metadata, check every node, envelope end,
composite budget and parent restriction consistency, and only then publish.
Composite budget uses the existing two-level application bound 5e-10 with
initial scale, final state and exchange scales; patch inventory and parent
consistency use its existing 2e-12 bound. No tolerances were changed.

The stream is caller-owned: position/close and scratch staging/file publication
are the caller's responsibility. A failed write can leave a partial stream;
this API does not replace a destination file. Reads consume only this envelope,
so a caller requiring EOF must check it. Optional selected identity/element
arguments pass through, but only the fixed-mechanism path was exercised here.
This is a restart record library, not yet a multilevel command-line application,
schedule identity envelope, filesystem-atomic writer or MPI restart path.

## Completed focused verification

GNU Fortran 13.3.0, `/usr/bin/gfortran`, existing Debug/Ninja
`build/selected-tree-mpi-debug`; no toolchain configuration change.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_eb_tree_checkpoint_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_eb_tree_checkpoint_3d$' --output-on-failure --output-log build/eb-tree-checkpoint-checks/ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **1/1 passed**, zero failures/skips, 0.38 s (total 0.42 s). A five-node,
three-level ratio-2 tree with planar z=0.31 geometry and nonuniform states
advances one HLLC step, round-trips through a scratch stream, and advances
again. Restored states/temperatures, every outer metadata field and the next
step's states/exchanges match uninterrupted execution exactly. Corrupting the
last envelope marker after valid patch records rejects without changing
any tree states/temperatures or the sentinel time. Wrong node count also
rejects transactionally. No selected-mechanism, filesystem publication or
changed-rank restart test was run.

Provenance, project contract and whitespace checks pass. No Python source
changed. No failed build/test attempt occurred. The test was rerun after
removing unused test declarations and improving failure messages.
Logs: `build/eb-tree-checkpoint-checks/{build,ctest}.log`.
