# EB particle field sampling and gas deposition — 2026-09-18

Implementation: `6361552345fb7f8c17d5ffb64f991923b40ae9ea`, based on `ce4bbcd`.

The single-grid API uses up to eight Cartesian cell-center weights, multiplies
them by fluid volume fractions, filters by open-face connectivity from the
containing cell within this support and normalizes. The original nodal level
set rejects solid query points. Sampling and deposition share these weights.
Integrated contributions are divided by actual fluid cell volumes, so their
volume integrals close. The gas depositor validates incoming active-cell states
and recovers final EOS temperatures before publishing both fields. Failed
queries or final states leave caller outputs unchanged; covered cells are not
updated or passed to the EOS.

## Completed focused checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; existing fixed tests-enabled and
selected H2/O2 tests-disabled configurations reused. Compiler
`/usr/bin/gfortran`, MPI wrapper `/usr/bin/mpifort`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_field --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_(field|drift|intersection)$' --output-on-failure --output-log build/spray-eb-field-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_mechanism_provenance.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **3/3 passed**, total **0.06 s**; field **0.01 s**, drift and intersection
each **0.00 s**. The new field fixture uses a nodal planar cut at x=1.25 in a
4-by-4-by-4 domain with unit cells. It checks normalized positive weights,
covered-cell exclusion, constant-field sampling, two integrated contributions
to N2 gas at 1000 K, every conserved component's fluid-volume inventory,
covered-state preservation, final EOS temperature, negative-density feedback
rollback and late solid-query sampling rollback. Inventory tolerance is 1e-14
relative to max(1, initial component inventory); temperature tolerance is 1e-8 K.
The fixture does not independently qualify arbitrary nonconvex connectivity.

Selected library compilation, provenance, Python syntax, project contract and
whitespace checks passed. No build/test failures or skips. Selected runtime
behavior was not rerun. Artifacts are in `build/spray-eb-field-checks/`:
`build.log`, `build-selected.log`, `ctest.log`.

Callers bind the nodal field to the supplied geometry. The stencil is a
volume-fraction-weighted approximation, not exact fluid-cloud intersection or
linear-exact interpolation near EB surfaces. The code has no small-cell
redistribution or source timestep control. Geometry validation is currently
repeated per query; production scaling is not qualified. Drag/heating/evaporation
stepping, distributed EB particles and EB application integration remain open.
No convergence, external solver parity or physical validation was performed.
