# Ambient molecular weight in vapor skin states — 2026-09-25

Implementation `7797e7ebcafb15dbab6d5adb54fe1e434da38a72`, base `95d75b9`.
Reference PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`, `Drag.H`
lines 127–162. SHA256
`bcc2bcb2f01b8cedfeb92d0bfc9ae2c8e723e77285f7011aeb77d5b138a29345` is pinned in
`references/spray_property_scope.json`; source bytes are retained under
`build/spray-property-reference-audit/`.

## Numerical correction

The reference vapor-mass denominator combines the liquid vapor mean weight
with the **complete** ambient gas mean weight. PeleF formerly used the carrier
mean weight excluding active fuel vapor. Both single-liquid and mixed-liquid
builders now use the complete gas composition. This matters in humid/fuel-bearing
gas and with independent liquid molecular weights.

The single-liquid builder now selects active vapor using the mole-pressure
criterion before constructing the skin state. It uses the liquid weight in
the vapor denominator and the mapped gas weight in the numerator. An active
surrogate can have surface mass fraction below ambient mass fraction; it still
participates in reference skin interpolation while the existing noncondensation
rate limiter prevents negative evaporation. Inactive vapor retains ambient
transport composition with zero active surface/ambient transfer fractions.
Single and mixed builders reject surface mass totals at or above one.

One-third droplet and one-half wall-film interpolation remain. Binary diffusivity,
heat-transfer coefficient conventions, conservative exchange and temporal
integration are unchanged; full reference coefficient/integrator parity is not
claimed. The regular context prefix is now `spray-v2-ambient-skin` and mixture
identity `mixed-liquid-components-v4-ambient-skin` reaches regular MPI, AMR and
EB context binding. Old numerical-model checkpoints are incompatible.

## Completed focused checks

GNU Fortran 13.3.0, Debug/Ninja, Python 3.12.3; unchanged fixed tests-enabled and
selected H2/O2 tests-disabled builds. Fixed core, regular, MPI AMR and MPI EB
applications built; selected regular/MPI AMR/MPI EB applications built.
Compiler warnings remain. No distributed execution was run for this increment.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core pelef_spray_3d pelef_mpi_amr_spray_3d pelef_mpi_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_3d_selected pelef_mpi_amr_spray_3d_selected pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_(vapor_skin|many_to_one|liquid_mw|mixture_rates|mixture_wall_rates|film_binary)$' --output-on-failure --output-log build/ambient-skin-checks/ctest.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/spray/test_spray_application.py --exe build/spray-amr-reactive-debug/pelef_spray_3d SprayApplication.test_molecular_weights_restart
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/spray/test_spray_application.py --exe build/spray-amr-reactive-selected-debug/pelef_spray_3d_selected SprayApplication.test_molecular_weights_restart
python3 -m py_compile tests/spray/test_spray_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Core **6/6 passed**, total **0.07 s**. The new test covers ambient water mass
fraction 0.2 with default/36/1000 kg/kmol liquid weights and an inactive case
at water mass fraction 0.9, at 373.15 K droplet temperature and 202650 Pa.
For both droplet and wall skin rules it checks reference-equation surface
fractions, activity, ideal-gas skin density and single/mixture equivalence at
2e-14 tolerance. Existing shared-vapor, independent-weight, mixture-rate,
wall-rate and binary-diffusion checks pass without tolerance changes.

Fixed/selected application methods each passed **1/1**, in **1.685/1.690 s**.
Each method covers scalar, port and mixed-liquid molecular-weight settings in
H2O/N2 ambient gas with mole ratio 1:9. Existing mass/energy/momentum checks,
exact gas/parcel CSV restart equality and changed/invalid-weight rejection pass.
These are 2-by-2-by-2 synthetic-liquid cases with gas chemistry disabled.

No failures or skips. Provenance, Python syntax, project contract and whitespace
passed. Artifacts: `build/ambient-skin-checks/{build-fixed.log,build-selected.log,
ctest.log,LastTest.log,application-fixed.log,application-selected.log}`.
No full matrix, installed-binary audit, external solver run, convergence or
experimental spray/film validation was performed.
