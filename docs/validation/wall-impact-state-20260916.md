# Splash children and deposited-film state — 2026-09-16

Source: parent `42bd4bd` plus this change. GNU Fortran 13.3 Debug,
`build/port-spray-debug`, MPI disabled.

## Scope

Construct four splash child parcels and a cylindrical stationary deposition
patch from supplied impact fractions. Child mass, velocity, placement and
multiplicity follow the pinned PelePhysics construction; film geometry follows
its modified Pasandideh-Fard expression. Child IDs/lineage/material identity and
reset deformation/RT state use PeleF conventions. A coordinate-axis fallback
resolves the reference tangent construction's doubly parallel case.

Separate signed wall impulse/energy budgets close each selected fraction's
inventory. Mass stays in children/film. Constructors leave the original parent
unchanged, so their caller must consume the selected fraction once. Invalid
inputs preserve destination state and return zero wall budget. This completes
state construction kernels, not CFD boundary events, parent transition dispatch,
film evolution, domain ownership or checkpoint integration.

## Focused evidence

Two CTest entries passed in 0.03 s, zero failures/skips: wall state and prior
wall correlations. New checks cover four equal-mass dry normal-impact children,
IDs/lineage/material/multiplicity, outward velocity and position, inventory plus
wall budget closure and invalid mass-fit rollback. Film checks cover deposited
mass/volume, wall momentum/energy, invalid Reynolds-number rollback and pinned
C++ geometry values.

For a 100-micrometer drop, 100 represented drops, deposited fraction .6,
density 1000 kg/m3, We=100, Re=1000 and contact angle pi/3, the isolated reference
C++ expression gives film diameter 0.0014595751300262862 m and thickness
1.8776173162247343e-5 m. Fortran agrees within 1e-13 relative.

Reference commit and hashes are in `references/spray_wall_impact.json`, including
`AhamedSplash.H` and newly recorded `SpraySB.cpp`. Build, mechanism provenance,
project contract and whitespace checks passed. No Python source changed and no
broad matrix or physical wall-impact validation was run.

Artifacts: `build/wall-impact-reference/film_oracle.{cpp,txt}` and its executable;
`build/port-spray-debug/wall-state-build.log`;
`build/port-spray-debug/Testing/Temporary/LastTest.log`.

```sh
c++ -std=c++17 build/wall-impact-reference/film_oracle.cpp -o build/wall-impact-reference/film_oracle
build/wall-impact-reference/film_oracle
cmake --build build/port-spray-debug --target test_spray_core --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_wall_(state|impact)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```
