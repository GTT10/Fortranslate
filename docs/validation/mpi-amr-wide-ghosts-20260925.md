# MPI AMR configurable halo width — 2026-09-25

Implementation `9e88bfa426fefe626f52c8e2b0c2fb019a61b5ed`, based on `cc279ba`.

`fill_mpi_amr_tree_ghosts_3d` now accepts `ghost_width`, default two. Effective
width enters collective min/max agreement alongside depth and endpoint presence;
different widths or values outside 1--128 reject before metadata walks. The
agreed width is forwarded through the shared sampler for counting, recording and
publishing donor requests, including coarse-time endpoint selection. Existing
owner-local output and atomic publication remain.

## Completed checks and runtime limitation

GNU Fortran 13.3.0 / Open MPI 4.1.6, Ninja Debug, existing configurations. Fixed
`test_mpi_amr_patch_tree_ghost_3d` and selected H2O2
`pelef_selected_mpi_amr_tree_reactive_3d_runtime` built successfully.
Provenance preflight, Python syntax, project contract and whitespace passed.

The test source adds seven-layer periodic serial/MPI comparison for PCM, MC and
minmod interpolation, owner-only output and rank-disagreeing-width rollback.
**These new MPI assertions have not been executed.** This is build evidence,
not numerical or MPI runtime qualification. The previously stalled CTest PID
2018710 and launcher PID 2018713 were still present at the initial observation;
no additional MPI execution was launched. Earlier serial wide-halo tests are
recorded separately and were not rerun here.

No runtime test successes, external parity, physical validation, scaling or full
configuration matrix are claimed for this increment. Dynamic MPI AMR LES
dispatch remains separate work.

## Commands and artifacts

Logs: `build/mpi-amr-wide-ghost-checks/` (`provenance.log`, `build.log`,
`build-selected.log`, `contract.log`).

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_mpi_amr_patch_tree_ghost_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_mpi_amr_tree_reactive_3d_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

After checking that the previous launcher is terminal, run the existing
`unit_mpi_amr_patch_tree_ghost_3d_2` focused test to verify the new assertions.
