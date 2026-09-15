# Regular-grid spray boundaries — 2026-09-15

Source: parent `130dbd3` plus this commit's spray modules, fixture, tests,
checker and docs. GNU Fortran 13.3 Debug, `build/port-spray-debug` with spray
enabled. This is focused algorithm/application evidence, not spray impact
validation or full PeleC equivalence.

## Implemented scope

Gas periodic/outflow/slip-wall/stationary adiabatic no-slip faces couple to
periodic wrapping, parcel outflow removal and ideal specular reflection.
Nonperiodic interpolation/deposition clamps at cell centers with eight distinct
indices and normalized weights. Reflection changes parcel normal momentum;
outflow removes its remaining inventory. Both join gas surface fluxes in the
transactional split-step ledger. Failed steps leave state, parcels and accepted
ledger unchanged. SGS must be zero for physical boundaries.

Checkpoint V2 stores the boundary ledger with the existing strict numeric
parser and atomic no-overwrite publication. The context binds all six faces.
Periodic checkpoints remain V1. Physical-boundary history adds five boundary
inventory columns; conservation residuals subtract them alongside injection.

## Focused checks

- `spray_boundaries`: physical stencil does not reference the opposite wall,
  reflection matches analytical normal impulse and conserves parcel energy,
  and outflow removes the complete parcel inventory.
- Existing `spray_integration_coupled`, `spray_integration_checkpoint` and
  `spray_integration_rollback` pass. Together with the boundary unit: 4/4.
- `spray_physical_boundary_restart`: wall and outflow application cases each
  compare uninterrupted runs with stop-after-two-steps/restart to 0.6 microseconds.
  Gas CSV, parcel CSV and final history row match exactly. All mass, momentum
  and energy residuals satisfy the checker's fixed bounds. Two wall parcels
  survive; both outflow parcels exit. Changing checkpoint face types is rejected
  before outputs are created.
- The final physical-boundary checker and existing `spray_application_cli`
  test passed 2/2 in 12.46 seconds.
- Project contract, pinned FFCM1_Red source provenance, Python syntax and
  `git diff --check` pass. No compiler/backend/MPI matrix was run.

The initial direct application artifacts are in
`build/port-spray-debug/boundary-check/spray-boundaries-qdkppy8o/`.
Wall mass/energy residuals there are `1.32135114931123467e-23` kg and
`-5.20995510474450998e-18` J; outflow residuals are
`1.12315617548215789e-23` kg and `-1.30104260698260532e-18` J.
CTest keeps subsequent artifacts beneath
`build/port-spray-debug/tests/spray_physical_boundaries/`.

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d test_spray_core test_spray_integration --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(boundaries|integration_(coupled|checkpoint|rollback))$' --output-on-failure
ctest --test-dir build/port-spray-debug -R '^spray_(physical_boundary_restart|application_cli)$' --output-on-failure
python3 -m py_compile tools/check_spray_boundaries.py
python3 tools/check_project_contract.py --source .
python3 tools/check_mechanism_provenance.py --bundle build/port-methanol-reference/ffcm1_red.json --source build/port-methanol-reference/ffcm1_red.yaml
```

Remaining: boundary-aware SGS, gas reservoirs, moving/isothermal spray walls,
sticking/splash/film models, multicomponent droplets, breakup/collision,
AMR/EB/MPI parcels and physical validation. Current integration cases use inert
water with drag; reacting/evaporating physical-boundary cases were not qualified
in this increment.
