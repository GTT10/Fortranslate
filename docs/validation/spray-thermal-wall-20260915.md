# Spray isothermal and moving gas walls — 2026-09-15

Source: parent `d55d9a8` plus this commit's wall plumbing, SGS flux, fixture,
tests and documentation. GNU Fortran 13.3 Debug, `build/port-spray-debug`.

## Implemented change

The spray application accepts prescribed wall temperature and tangential wall
velocity, passing them through the coupled operator to molecular transport and
SGS. Molecular heat/work reuses the existing regular-grid implementation.
SGS ghosts reflect velocity around the prescribed wall velocity. SGS conduction
uses the wall-to-cell temperature difference over half a cell, with doubled
boundary diffusion allowance in its stability check. All face fluxes contribute
to the same boundary ledger; failed steps retain transactional rollback.
Isothermal/moving application walls require molecular transport. Normal wall
motion rejects. Active wall parameters bind the checkpoint context.

Particle walls remain specular. This is fixed-boundary gas heat/momentum
exchange, not a splash, adhesion, particle-wall heat transfer or moving geometry
model. SGS still uses the constant-coefficient closure without wall damping.

## Focused evidence

Six tests passed in 16.00 seconds:
`spray_thermal_wall_0`, `spray_thermal_wall_.17`, `spray_les_boundary`,
`spray_physical_boundary_restart`, `spray_sgs_boundary_restart`,
`spray_application_cli`.

The SGS unit checks wall mechanical work equals wall speed times traction
within relative 1e-12. A hotter wall adds heat and its energy ledger matches
volume-integrated divergence within 1e-12. Normal wall motion rejects without
changing gas state or temperature.

The coupled fixture uses a 4x4x4 grid, nitrogen at 900 K and 101325 Pa,
1000 K walls, a 10 m/s tangential top wall, and two reflecting water parcels.
At 0.6 microseconds, both parcels remain live. With SGS off, net boundary energy
is 7.3575454248e-7 J and energy residual -6.2596e-18 J. With Cs=.17 and imposed
initial shear, these are 3.9914461512e-6 J and 4.4800e-18 J. These two fixtures
have different initial shear and do not isolate the effect of Cs.

Stopping after two steps and resuming gives exact gas/parcel CSV and final
history agreement in both cases. Changing wall speed from 10 to 11 m/s rejects
the checkpoint before creating output files. Existing wall/outflow and CLI
checks also pass. This does not establish cross-commit bitwise parity.

Artifacts under `build/port-spray-debug/tests/`:
- `spray_thermal_wall_0/spray-boundaries-mehn5ps9/`
- `spray_thermal_wall_.17/spray-boundaries-zkrdueir/`

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d test_spray_core test_spray_integration --parallel 8
ctest --test-dir build/port-spray-debug --output-on-failure -R '^spray_(thermal_wall_.*|les_boundary|physical_boundary_restart|sgs_boundary_restart|application_cli)$' -j 4
python3 tools/check_mechanism_provenance.py --bundle build/port-methanol-reference/ffcm1_red.json --source build/port-methanol-reference/ffcm1_red.yaml
python3 tools/check_project_contract.py --source .
python3 -m py_compile tools/check_spray_boundaries.py
git diff --check
```

Provenance, project contract, Python syntax and whitespace checks pass. Initial
builds caught a private shared validator and an overlong namelist line; both
were corrected. No compiler/backend matrix or physical wall/LES validation was run.
