# Serial physical-boundary AMR regridding — 2026-09-17

Source: parent `0396b7c` plus this change. Logs and source hashes:
`build/amr-boundary-regrid-checks/`. Reused GNU Fortran 13.3.0 Debug build.

The conservative 3D remap API now accepts physical boundary data. It restricts
the old hierarchy, fills physical coarse ghosts and uses those neighbors for
limited conservative prolongation at newly refined boundary cells. Positivity
backtracking, species closure, exact overlap reuse, final average-down and
transactional commit retain their existing behavior. Periodic seam patches are
still rejected. The temperature tagger can include physical boundary cells and
clip its buffer to those faces; periodic axes retain an interior-cell margin.

Serial applications opt into this behavior with
`regrid_physical_boundaries = .true.`. The default preserves the former interior
tagging/remap path, including existing serial/MPI dynamic cases. MPI frontends
reject the new option until their boundary-touching stepping is connected.
Dynamic checkpoints include `PHYSICAL_TAGGING_V1` in their existing signature
record when enabled. The old signature remains unchanged when disabled.
Restored geometry may touch physical faces only with the matching tagging mode;
periodic seam geometry remains invalid. No checkpoint schema number changed.

The remap unit checks physical corner expansion, shrinkage and movement,
conserved composite inventories, exact overlap, inclusion of boundary tags,
invalid-reservoir rollback and periodic-seam rejection. Existing interior linear
prolongation/tagging tests remain. Application regressions expand wall/corner
patches at step zero, checkpoint, and resume full-H2/O2 and elementary reacting
cases. Coarse/fine CSVs and boundary-ledger output match uninterrupted execution
exactly. The existing context mismatch probes remain in the regression.

Commands:

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
cmake --build build/spray-mpi-migration-gnu-debug --target pelef_amr_reactive_3d test_amr_regrid_3d test_amr_reactive_3d_config test_amr_reactive_3d_checkpoint --parallel 4
ctest --test-dir build/spray-mpi-migration-gnu-debug -R '^(unit_amr_regrid_3d|unit_amr_reactive_3d_(config|checkpoint)|amr_boundary_touching_dynamic_restart_3d)$' --output-on-failure
python3 -m py_compile tools/check_amr_physical_boundary_restart.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Result: 4/4 tests passed, no skips (config 0.00 s, checkpoint 0.09 s,
remap 0.15 s, application 53.44 s; total 53.73 s). Provenance, Python syntax,
project contract and diff checks passed. The preliminary run before adding
the explicit tagging-mode setting is kept separately and is not counted here.
Final application artifacts:
`build/spray-mpi-migration-gnu-debug/tests/amr_boundary_touching_dynamic_restart_3d/run-5dnfpfgx/`.
Runnable example: `cases/amr_reactive_3d/dynamic_touching_wall.nml`.

Boundary-touching MPI stepping/remapping and general multipatch/multilevel
hierarchies remain open. This is numerical conservation/restart evidence,
not physical validation.
