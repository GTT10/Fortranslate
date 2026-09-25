# Fine-resolution EB geometry halos — 2026-09-25

Code `cf9327a`, parent `1e12143`. Extracted EB patches now retain one cell
of same-resolution volume fractions, including edge/corner context. Extraction
copies available source cells and propagates existing context through nested
extraction. Missing exterior context is marked -1, not invented by clamping.
Geometry validation checks bounds, finite admissible values and interior
agreement. Builders and coarsening still produce ordinary geometry without
this optional context.

Two-level/multipatch parent halo fill, serial common/recursive tree halo fill
and MPI tree halo fill prefer retained target-resolution fractions. State and
temperature donors retain their existing same-level/ancestor and temporal
interpolation policy. Physical-domain coordinates are clamped before looking
up metrics. Covered fine ghosts do not require EOS recovery. Geometry created
without exterior context retains the historical donor-fraction fallback.
This removes the parent-constant fraction approximation for the standard tree
initialization/remapping paths that extract from full level geometry; it does
not claim arbitrary missing geometry can be reconstructed from a coarse cell.

## Restart contract

The generic EB checkpoint writes an optional `GEOMETRY_HALO_V1` section before
`DIAGNOSTICS`. Context is restored and compared against the expected geometry
before publishing state. Records without this section retain their existing
layout and remain readable when expected geometry also has no context. A legacy
cropped-patch record without context is rejected when the current expected
patch has context; silently discarding or inventing it would change continuation.
Tree wrappers inherit this record extension. Changed-rank restart was not rerun.

## Completed checks

GNU Fortran 13.3.0, OpenMPI 4.1.6, Python 3.12.3, reused Ninja Debug build
`build/spray-amr-reactive-debug`. Four distinct CTests passed: standalone EB
checkpoint 0.14 s, positive-viscosity two-level 1.12 s, tree checkpoint 0.38 s,
three-rank positive-viscosity tree 1.85 s (3.54 s total). After adding the
context-mismatch assertion, the tree checkpoint alone passed again in 0.38 s
(0.43 s total). No runtime failures or skips.

The two-level test compares every exterior fraction directly to the full fine
geometry and requires a value differing from its coarse donor, for ratios 2/3
and physical/interior patches. Existing common/recursive conservation and
rollback checks remain. MPI common hydro retains serial state/temperature and
ledger equivalence under distributed and single-owner layouts. Tree restart
checks exact context, state, temperature and subsequent update/ledger equality,
plus context mismatch and late-read rollback. No tolerance or golden changed.

Two initial compilation failures were corrected: a geometry lookup was inserted
in the peer routine rather than parent routine, and a test expression had an
extra parenthesis. Their logs are retained. Final builds, mechanism provenance,
Python syntax, project contract and whitespace checks passed. Compiler warnings
remain. Full configuration matrix, recursive MPI runtime, external EB geometry
parity and physical validation were not run.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_two_level_3d test_mpi_eb_tree_hydro_3d test_eb_tree_checkpoint_3d test_reactive_eb_3d_checkpoint --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^(unit_eb_two_level_viscosity_3d|unit_mpi_eb_tree_viscosity_3d_3|unit_eb_tree_checkpoint_3d|unit_reactive_eb_3d_checkpoint)$' --output-on-failure
# Final checkpoint-only rebuild/recheck after adding context mismatch coverage:
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_eb_tree_checkpoint_3d --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^unit_eb_tree_checkpoint_3d$' --output-on-failure
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Artifacts: `build/eb-geometry-halo-checks/` contains `build.log`, `rebuild.log`,
`final-build.log`, `test-build.log`, `checkpoint-build.log`,
`LastTest-first.log` and `LastTest-checkpoint.log`.
