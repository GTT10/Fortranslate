# Multilevel EB automatic common intervals — 2026-09-18

Implementation: `7fc10832853c5d95e385400f23095c54637f28c9`, based on
`806019a`. Library API: `reactive_eb_tree_timestep_3d`.

The estimator averages a private snapshot, gathers maximum Cartesian/wall
signal speeds and molecular diffusivity across all fluid-containing nodes,
and passes these donor bounds to existing per-patch Euler hydro and diffusion
interval routines. It takes the minimum across all nodes, including hidden
parent cells advanced before restriction. Empty-fluid nodes do not constrain
the estimate. Source tree and output interval remain unchanged on failure.
The full diffusion bound is conservative for the coupled half steps; no
chemistry accuracy estimate or redistribution is introduced.

The application adds `automatic_timestep`, `hydro_cfl`, `transport_cfl`, caps
the result by `time_step` and remaining time, and retains bounded retry halving.
V2 application envelopes store interval mode and CFL coefficients. V1 parsing
remains restricted to fixed mode with existing embedded configuration checks;
legacy V1 compatibility was not rerun here. V2 policy changes reject restart.

## Completed focused verification

GNU Fortran 13.3.0, `/usr/bin/gfortran`, existing Debug/Ninja
`build/selected-tree-mpi-debug`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tools/check_eb_tree_application_3d.py
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target pelef_reactive_eb_tree_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^regression_eb_tree_(application|auto_interval)_3d$' --output-on-failure --output-log build/eb-tree-auto-checks/ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

Final **2/2 passed**, zero failures/skips, total 18.28 s (fixed 2.45 s,
automatic 15.80 s). The automatic fixture uses hydro CFL 1e-5 to activate the
bound below its 1e-9 cap: first interval **1.3634407285670431e-10**, 15 accepted
steps to 2e-9, minimum interval **9.1183116651589353e-11** including the final
clip. The maximum normalized application inventory residual is
**3.4575864181991328e-15**. Both modes preserve exact uninterrupted/restarted
state and budget CSVs, independent CSV volume-integral checks, and protected
output files on rejected policy/mode or corrupted-record restarts.

Initial build rejected an overlong namelist line; continuation fixed it.
An accidentally started check then exercised the stale executable, which
rejected the new input. The first current-binary automatic fixture at CFL
1e-4 was capped at 1e-9 and failed the test's requirement to exercise actual
shortening. Only that fixture CFL was reduced to 1e-5; no solver gate changed.
Those attempts are retained as `initial-build.log`, `stale-binary-ctest.log`
and `cap-limited-ctest.log` in `build/eb-tree-auto-checks/`. Final build/test
logs are in the same directory. Runtime artifacts are under
`build/selected-tree-mpi-debug/tests/eb_tree_auto_interval_3d/`.

Provenance, Python syntax, project contract and whitespace checks pass.
This is a bounded application/restart check, not temporal convergence,
large-tree performance, transport-dominated stability or physical validation.
