# Coupled uniform acceleration — 2026-09-15

Source: parent `7435354` plus this commit's uniform-force operator, application,
source ledger/checkpoint, tests and documentation. GNU Fortran 13.3 Debug in
`build/port-spray-debug`.

## Implementation

A constant user acceleration acts on both gas and liquid parcels. Exact
velocity kicks surround the existing coupled step. Gas energy receives
`(momentum_old + delta_momentum/2) dot delta_velocity`; parcel work uses its
represented mass. These increments accumulate separately from interphase
feedback and boundary flux. Gas temperature is preserved and thermal recovery
checked. Failed stages leave the original gas, parcels and source increments
unchanged. An application timestep cap limits acceleration displacement.

The active acceleration vector binds restart context. `_FORCE_V1` extends each
existing checkpoint format with the external source ledger; default zero-force
formats remain unchanged. History reports external impulse/work and subtracts
it in conservation residuals. No potential energy, hydrostatic well-balancing,
spatial gravity or phase-specific forcing model is implied.

## Focused evidence

Four distinct tests passed: `spray_integration_acceleration`,
`spray_acceleration_restart`, `spray_integration_checkpoint`,
`spray_application_cli`.

The periodic force unit disables drag, heat, evaporation, SGS and chemistry.
At acceleration (0,0,-9810) m/s² over 0.2 microseconds, particle velocities
match v0+a*dt within absolute 1e-13 m/s, positions match x0+v0*dt+a*dt²/2 within
1e-16 m. Gas velocity matches the same acceleration; temperature changes by
less than 1e-8 K. External momentum balance is within absolute 1e-18 kg m/s
and energy within 1e-12 of initial energy. A subsequent invalid chemistry
backend rejects after a force kick, preserving both phases and zero returned
source increments.

The wall/outflow CLI fixtures retain drag and use acceleration -9810 m/s²,
stop after two steps, and resume to 0.6 microseconds. Gas/parcel CSVs and final
history rows match uninterrupted runs exactly. Changing acceleration rejects
restart before output creation. The wall run external z impulse is
-1.4348302131e-10 kg m/s, external work 5.1048926117e-11 J, and energy residual
-4.0905e-18 J. An artifact set is
`build/port-spray-debug/tests/spray_acceleration/spray-boundaries-0s7ljihw/`;
the final rerun after adding thermal recovery uses another unique directory.

The first application restart check failed because its ledger consistency
check omitted external work. Adding that source term fixed it; no tolerance
was relaxed. Initial compilation also caught an overlong signature and it was
wrapped. The final force unit/restart pair passed after thermal recovery was
added. Existing zero-force checkpoint and CLI checks passed separately.

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d test_spray_integration --parallel 8
ctest --test-dir build/port-spray-debug --output-on-failure -R '^spray_(acceleration_restart|integration_(acceleration|checkpoint)|application_cli)$' -j 4
python3 tools/check_project_contract.py --source .
python3 tools/check_mechanism_provenance.py --bundle build/port-methanol-reference/ffcm1_red.json --source build/port-methanol-reference/ffcm1_red.yaml
python3 -m py_compile tools/check_spray_boundaries.py
git diff --check
```

Project contract, FFCM1_Red provenance, Python syntax and whitespace pass.
No compiler/backend matrix, hydrostatic benchmark or physical settling
validation was performed.
