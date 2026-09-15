# Finite circular nozzle outlets — 2026-09-16

## Implementation and provenance

Tested parent: `13ec0807ee5e40116eec2f844486fbced41b616a` plus this commit's source/tests. Exact tested
hashes and logs are in `build/nozzle-checks/`.

`references/spray_nozzle.json` pins PelePhysics
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`, `SprayInjection.H` and `SprayJet.H`,
with their hashes and license. A positive `nozzle_radius` now distributes
parcel positions over a circular outlet plane. The source's linear radius
sampling is retained with deterministic midpoint quantiles and ID-based
azimuth. Optional per-port radii connect to initial and shared/independent
scheduled injection. Zero keeps the existing point-source calculation.

The application validates the complete disk's Cartesian extents, rather than
only its center or sampled points, before running. Radius changes enter
checkpoint identity. Particle quantities and their conservative accounting
are retained; this is not a new gas boundary condition.

## Focused verification

GNU Fortran 13.3.0, existing Debug spray build, MPI disabled:

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d test_spray_core --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(injection|application_cli)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Two CTest entries passed in 36.51 s, including 28 CLI tests; no failures or
skips. Build, mechanism provenance, Python syntax and project contract passed.
Compiler warnings include real comparisons and existing bounds-analysis
messages.

The new cases check a rotated single outlet, two different outlet radii/axes,
and independently scheduled outlets. With stationary injected parcels, all
positions satisfy the plane equation and expected linear radial quantiles.
Cases also check system conservation, byte-identical gas/parcel CSVs and
identical final history after split restart. Changed radius rejects restart;
negative radius, a disk outside the domain and unused port settings fail.
Existing point-injection checks also pass.

## Remaining boundary

The sampling is uniform in radius, not area. Deterministic joint quadrature
is not an independent-random physical trajectory model. Hollow jets, swirl,
time-varying outlet geometry, gas inlet qualification and physical spray
validation remain open. No full-suite or external solver run was performed.
