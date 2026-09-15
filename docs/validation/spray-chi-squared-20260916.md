# ChiSquared diameter table — 2026-09-16

## Source and implementation

Tested parent: `f610fdbe8d7e1d5e1e10650331a15d58769f40ae` plus this commit's
source/tests. Tested hashes and logs are retained in `build/chi-size-checks/`.

`references/spray_size_distribution.json` pins PelePhysics
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`,
`Source/Spray/Distribution/Distributions.cpp`, and its license.
`diameter_distribution='chi-squared'` now implements the source's 100-bin
cumulative table and left-edge diameters. `drop_diameter` maps to nominal d32,
including per-port overrides. The table is conditioned on the configured
positive diameter bounds, removing the source's zero-diameter first bin.
No supported bins is an error. Midpoint quantiles, mass/number weighting and
exact-dose normalization reuse the injection path. The existing selector
context distinguishes checkpoints from other distribution laws.

## Focused checks

GNU Fortran 13.3.0, existing Debug spray build, MPI disabled:

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d test_spray_core --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(size_distribution|application_cli)$' --output-on-failure
ctest --test-dir build/port-spray-debug -R '^spray_application_cli$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```

The first run passed the size-distribution kernel test but failed the new CLI
exact upper-bound assertion: CSV reconstructs diameter from mass and reported
2.00000000000000152e-5 m for a 2e-5 m boundary. The comparison now permits
4e-15 relative roundoff; no model or other test tolerances changed.
The subsequent CLI run passed all 26 tests, with no skips. Build, mechanism
provenance, syntax and contract checks passed. Compiler warnings include real
comparisons and existing bounds-analysis diagnostics. Initial and final test
logs are retained separately.

The new case checks the eight positive-support midpoint bins
[12,18,23,28,33,38,46,61] for two different nominal d32 values, equal represented
counts, total conservation, byte-identical gas/parcel CSVs and identical final
history after restart. It checks explicit diameter truncation, empty support,
unused standard deviation and changed-model restart rejection. Shared-port
routing is exercised directly; independent schedules reuse the same injector.

## Boundary

The pinned discrete law is retained; a continuous chi-squared inverse is not
substituted. Deterministic quantiles, positive-support conditioning and finite
mass normalization differ from the reference random sampler. Nominal d32 is
not a guarantee of the measured post-truncation Sauter mean. Physical child
creation/dispersion, mixture droplets and particle MPI/AMR/EB remain open.
No full-suite or external solver/physical validation was performed.
