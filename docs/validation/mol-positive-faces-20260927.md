# MOL positive-face slope bound — 2026-09-27

The [particle diagnostic](mpi-mol-particle-diagnostic-20260927.md) identified
negative trace-fuel partial density from separately limited characteristic
amplitudes. `limit_mol_positive_slope_3d` now scales all amplitudes of one cell
by the same factor, bounded by both reconstructed partial densities and
pressure. It verifies both actual floating-point faces before publishing the
slope. An active bound is moved slightly inward; any remaining rounding
failure reduces the factor again. Invalid cell centers still reject, and the
raw face reconstruction still rejects negative species. Cell averages are
untouched; both neighbors continue to share one conservative interface flux.

This is an additional positivity safeguard beyond pinned PeleC MOL, not a
claim of exact reference reconstruction for affected faces. Smooth admissible
slopes remain unchanged. No species acceptance tolerance or conservation
threshold was relaxed. The diagnostic's very small negative fuel density is
not accepted or clipped to zero.

## Completed checks

GNU Fortran 13.3.0 / OpenMPI 4.1.6, existing selected heptane Simple-transport
Debug MPI/spray build, tests/SUNDIALS disabled:

```sh
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2 > build/mol-positive-build.log 2>&1
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/heptane-simple-mpi-debug --target pelef_core --parallel 2 > build/mol-positive-core-build.log 2>&1
env -u LD_LIBRARY_PATH -u LIBRARY_PATH /usr/bin/gfortran -O0 -g -fcheck=all -ffpe-trap=invalid,zero,overflow -ffree-line-length-none -I build/heptane-simple-mpi-debug/modules tests/unit/test_mol_reconstruction_3d.F90 build/heptane-simple-mpi-debug/libpelef_core.a -o build/mol-positive-unit/check
build/mol-positive-unit/check > build/mol-positive-unit/run-final.log 2>&1
python3 build/mpi-mol-particle/run_positive_final.py > build/mpi-mol-particle/positive-final.log 2>&1
python3 -m py_compile build/mpi-mol-particle/run_positive_final.py
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```

These checks passed. The reconstruction executable passes 47 checks including
the exact previously failing center/slopes, both faces in three directions,
nonzero bounded slopes, unchanged smooth reconstruction, pressure bounds and
transactional invalid-center rejection. The first compile attempted to use
the stale fixed-core module after only the selected application build and
failed; rebuilding `pelef_core` resolved that build dependency.

The final application driver completes seven launches: common MOL continuous,
checkpoint, one-to-two-rank restart, MOL-to-CTU restart rejection, explicit
recursive-particle-MOL rejection, CTU checkpoint, and CTU-to-MOL rejection.
The negative cases must exit unsuccessfully with the designated message and
no final CSV. The original 4³ heated-wall fixture with four evaporating parcels
now completes two `1e-7 s` steps. Gas CSV bytes and sorted particle rows match
across restart. The largest absolute combined budget residual is
`2.710505431213761e-19` (components have different units; this is not a relative
physical-error measure). Detailed commands and budgets are recorded in
`build/mpi-mol-particle-positive-final/result.json`; binary SHA256 is
`b1300c674a331cb2b710fc57acb932bfae6bd46dd7827a9d7287c7960d89b20b`.

An earlier driver failed its recursive-rejection expectation because a later
`subcycling=.false.` namelist entry overrode its inserted true entry. The
corrected final driver fixes that duplicate value and uses a fresh directory.
Earlier artifacts remain in `build/mpi-mol-particle-positive`. Nonfatal
array-temporary, X11 authorization and executable-stack warnings remain.

Recursive particle MOL, reacting particle comparison, full-duration Spray-EB,
full configuration matrices and upstream temporal parity were not completed.
The separate CTU Jet failure is not repaired by this MOL change.
