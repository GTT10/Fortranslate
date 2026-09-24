# Watson latent heat and PeleMP boiling kernels — 2026-09-25

Implementation `e70a3b55a091f39103441d883e1b05f08b0eb95e`, base `7f83398`.
Reference PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`;
the five source hashes and call-site scope are recorded in
[the property audit](spray-property-scope-20260925.md) and
[`spray_property_scope.json`](../../references/spray_property_scope.json).

## Implemented boundary

`spray_parcel_mod` now exports two pure scalar property APIs.
`watson_latent_heat` evaluates the pinned correlation for a supplied latent
reference state and a subcritical evaluation temperature.
`pelemp_boiling_temperature` evaluates that correlation at a separately supplied
normal boiling point, then applies the pressure-dependent Clausius inverse
temperature and critical-temperature cap from `SprayData::calcBoilT`.

The normal boiling reference is 101325 Pa. Temperatures use K, latent heat J/kg,
pressure Pa and molecular weight kg/kmol. As in existing property guards,
temperatures are bounded by 1–100000 K, reference latent heat by 1e10 J/kg,
pressure by 1e-6–1e12 Pa and molecular weight by 1e-3–1e5 kg/kmol.
The latent reference and normal boiling temperatures must be below Tc.
Invalid/nonfinite inputs, an overflowing division, nonpositive inverse boiling
temperature or a resulting temperature below 1 K return zero with `ok=.false.`.
This deliberately rejects the reference formula's nonphysical negative branch.

No liquid record, application input, stored liquid energy or restart format
changed. The existing constant-latent surface cap is still the runtime path.
Connecting these APIs to that path requires distinct inputs and restart binding;
ordinary EOS vapor-state latent heat remains a separate enthalpy calculation.

## Focused checks

GNU Fortran 13.3.0, Debug/Ninja; unchanged fixed tests-enabled and selected H2/O2
tests-disabled configurations. Python 3.12.3 supplied a 60-digit Decimal oracle.
The new test compares four temperatures and four pressures with independent
decimal constants at relative tolerance 2e-14. It deliberately uses Tref=300 K
and normal Tb=373.15 K, with Lref=2.5e6 J/kg and Tc=647 K.
The pressure cases span 0.1, 1, 10 and 1000 atmospheres, including the Tc cap.
The test also checks the exact latent reference identity, NaN in each input,
invalid reference/critical states, zero latent/pressure, negative inverse
temperature and a tiny denominator under floating-point traps.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_runtime --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_(watson_properties|surface_cap|evaporation)$' --output-on-failure --output-log build/watson-properties-checks/ctest.log
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Both builds passed; selected runtime was compiled, not executed. CTest **3/3
passed**, total **0.05 s**, with no failures or skips. Existing surface-cap and
evaporation tests passed without tolerance changes. Provenance, Python syntax,
project contract and whitespace checks passed. Build warnings remain (277 in
each build log); this is not a warning-free build claim.

Artifacts: `build/watson-properties-checks/{build-fixed.log,build-selected.log,
ctest.log,LastTest.log,decimal_oracle.py,decimal_oracle.txt}`. The oracle script
records how the independent constants were produced. No new coupled application
execution was warranted because the new API is not yet dispatched there.
No full matrix, external solver run, install audit, convergence or physical
spray validation was performed. The tests establish the isolated formulas and
guards, not full PeleC spray parity.
