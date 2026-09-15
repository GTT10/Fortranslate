# Nozzle profiles and swirl — 2026-09-16

## Source and implementation

Tested parent: `c7ac0debcc8c858b565a48ed6570614a5af8ffe2` plus this commit's source/tests.
Exact tested file hashes and build/test logs are in `build/nozzle-profile-checks/`.

`references/spray_nozzle.json` pins PelePhysics SprayInjection.H, SprayJet.H
and SprayJet.cpp at `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b` with hashes
and license. The `solid` profile connects the source's radial inclination
relation; `hollow` connects rim placement and inclination spread. Swirl rotates
velocity azimuth relative to position azimuth. The existing cone remains the
`legacy` default. Positive radius is required for solid/hollow.

Global/per-port records route through initial, shared and independent scheduled
injection and checkpoint identity. Optional complete `spray_nozzle_model`
records expose the same path to Fortran callers. Angles use radians, with
explicit per-field inheritance and unused-port rejection in the CLI.

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

Two CTest entries passed in 38.86 s, including 29 CLI tests; no failed or
skipped tests. An initial build encountered the Fortran 132-column limit in
an extended argument list; line wrapping fixed it before the successful build.
Compiler warnings include real comparisons and existing bounds diagnostics.
Provenance, syntax and project contract passed.

The new checks recover the initial outlet location from ballistic motion and
verify rim/linear radii, axial speed, total speed, and radial/tangential
velocity components against the prescribed spread/swirl equations. Single,
shared two-port and independently scheduled two-port cases check conservation
and byte-identical final gas/parcel CSVs after restart. Changed swirl rejects
restart; zero radius for hollow, excessive swirl and unused port models fail.
Legacy nozzle/injection tests remain covered.

## Boundary

Deterministic midpoint spread/azimuth quadrature replaces the source's random
sampling; this is not stochastic trajectory parity. Wide hollow spreads can
produce negative or backward inclinations, as in the reference formula. Gas
boundary conditions remain separately configured. Time-dependent nozzle
geometry, physical breakup-child creation/dispersion, mixture droplets and
particle MPI/AMR/EB remain open. No full-suite or external physical validation
was run.
