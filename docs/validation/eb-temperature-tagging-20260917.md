# Exposed-fluid EB temperature tagging — 2026-09-17

Implementation: `d75b1f2270319933dfee3e7722b124ef832eef76`, based on
`fca69ed`. Adds `tag_reactive_eb_temperature_3d` to the two-level EB library;
existing update and remap algorithms are unchanged.

## Contract

The API receives coarse geometry/temperature, old fine geometries and
coarse-index origins, refinement ratio, positive finite temperature threshold,
coarse tile dimensions, nonnegative buffer width and maximum patch count.
It returns allocatable 3-by-count inclusive coarse-index lower/upper bounds.

Only exposed fluid cells vote. Fine hot cells are mapped directly to their
coarse parent indices, without averaging away fine temperature peaks. Coarse
cells hidden by fine patches and solid cells are ignored even if their stored
temperatures are stale or nonfinite. Every exposed fluid temperature must be
finite and positive. Fine geometries must coarsen to matching parent metrics;
overlap, malformed dimensions or invalid geometry reject.

The inclusive threshold and deterministic k/j/i tile ordering follow the
regular AMR tagging convention. A clipped coarse-cell buffer surrounds tagged
cells. Tiles partition the domain without overlap, include partial edge tiles,
and omit all-solid tiles. No hot cells returns zero patches. Exceeding capacity
rejects rather than silently dropping hot regions. Output bounds publish only
after complete success; failures preserve prior output arrays. No periodic
wrapping or additional threshold hysteresis is introduced.

## Completed focused checks

GNU Fortran 13.3.0, `/usr/bin/gfortran`, Debug/Ninja, existing
`build/selected-tree-mpi-debug`.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/selected-tree-mpi-debug --target test_eb_temperature_tagging_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/selected-tree-mpi-debug -R '^unit_eb_temperature_tagging_3d$' --output-on-failure --output-log build/eb-temperature-tagging-checks/ctest.log
python3 tools/check_project_contract.py --source .
git diff --check
```

**1/1 passed**, 0.04 s (total 0.08 s), zero failures/skips. A 4-cubed coarse
mesh, planar z wall at 0.31, and ratios 2/3 exercise a fine-only hot cell at the
exact threshold; ignored stale hot/nonfinite covered parents; patch permutation;
no-hot zero-patch output; a domain-corner buffered tag with explicit expected
bounds; capacity failure and rollback; invalid exposed temperature and rollback;
an integer-maximum buffer clipped safely to the domain; all-solid exclusion;
and partial edge tiles with explicit expected bounds.

Provenance, contract and whitespace checks passed. No Python files changed.
Logs: `build/eb-temperature-tagging-checks/build.log` and `ctest.log`.

This completes the tagging library operation, not adaptive application
execution. The application still needs temperature-policy inputs, zero-patch
handling, tag/remap dispatch and adaptation-bound restart. General multilevel,
MPI, nonplanar tagging qualification, external parity and physical/convergence
validation remain open.
