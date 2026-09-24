# Independent liquid molecular weights — 2026-09-25

Implementation `b2b8908c8bb1a7a8231fb9af03314709e61f66d1`, base `aa2c8e3`.
Reference: PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`.
`Source/Spray/SprayProperties.H` was extracted with `git show` at that revision;
SHA256 `9fc95260ec730e91d9e4555ad3000e5084bc5e31edf2329f1efa180fabc22550` is added
to `references/spray_property_scope.json`. Bytes are retained in
`build/liquid-molecular-weight-checks/SprayProperties.H`.

## Reference and implementation

`SprayProperties.H:71–147` uses gas molecular weights for one-to-one mapping,
reads user liquid molecular weights for many-to-one mapping and converts mass
to mole fractions using those liquid values. The reference labels many-to-one
mapping as under development. `Drag.H:127–142` uses the liquid vapor mean
weight in its denominator and mapped gas-species weights in vapor numerators.

PeleF liquid records now accept a molecular-weight override in kg/kmol. Zero
retains the gas-weight fallback. Nonzero weights are finite and bounded by
1e-3–1e5 kg/kmol. The resolved value feeds Raoult liquid mole fractions,
Clausius saturation pressure, legacy/PeleMP boiling caps and wall-impact boiling.
Mapped gas molecular weights still form the surface vapor mass numerator.
Gas thermodynamics and gas/parcel energy conventions are unchanged.

Regular scalar, port and component settings and both single-level/patch-tree
EB and MPI AMR component namelists bind the new values to restart context.
Temporary mixture records retain valid weights; their mechanical use does not
alter component inventory. The mapped gas species continues to define species
and elemental bookkeeping. Arbitrary real-liquid elemental formulas are not
introduced by a molecular-weight override.

The impact helper also now dispatches `pelemp-cap` to the Watson-based boiling
kernel. Previously only saturation used that branch. Regular PeleMP contexts
now include an impact-model marker; expanded MPI AMR/EB namelists also reject
older contexts. Cross-build restart compatibility is not promised. Legacy
regular inputs with zero molecular weights and without PeleMP keep their
existing context representation.

## Focused evidence

GNU Fortran 13.3.0, Debug/Ninja, Python 3.12.3. Existing fixed tests-enabled and
selected H2/O2 tests-disabled builds were reused. Fixed core, regular spray,
MPI AMR spray and MPI EB spray built; selected versions of all three applications
built. MPI applications were not executed for this change. Compiler warnings
remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core pelef_spray_3d pelef_mpi_amr_spray_3d pelef_mpi_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_3d_selected pelef_mpi_amr_spray_3d_selected pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_(liquid_mw|many_to_one|pelemp_cap|mixture_wall_event)$' --output-on-failure --output-log build/liquid-molecular-weight-checks/ctest-final.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/spray/test_spray_application.py --exe build/spray-amr-reactive-debug/pelef_spray_3d SprayApplication.test_molecular_weights_restart SprayApplication.test_pelemp_cap_restart
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/spray/test_spray_application.py --exe build/spray-amr-reactive-selected-debug/pelef_spray_3d_selected SprayApplication.test_molecular_weights_restart
python3 -m py_compile tests/spray/test_spray_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final core **4/4 passed**, total **0.06 s**. Equal liquid mass fractions with
weights 18/36 yield mole fractions 2/3 and 1/3; at the reference temperature
and two reference atmospheres the surface values are 1/3 and 1/6. The new test
checks those values, the Clausius exponent, mapped gas-weight numerator,
diffusion weights, impact boiling and invalid weights. Formula tolerances are
2e-14; existing default shared-vapor, PeleMP cap and wall-event tests also pass.

Fixed application methods **2/2 passed in 3.431 s**; selected **1/1 passed in
1.711 s**. Molecular-weight checks cover scalar, per-port and shared-vapor
mixture inputs, existing mass/energy/momentum bounds, exact gas/parcel CSV
restart equality, changed-weight restart rejection and invalid-weight rejection.
The fixed run also repeats the three PeleMP restart variants after the new
impact-context marker. Fixtures use synthetic liquid inputs with chemistry off.
No failures, skips or loosened tolerances. Provenance, Python syntax, project
contract, whitespace and manifest JSON parsing passed.

Artifacts: `build/liquid-molecular-weight-checks/`, including final fixed/selected
build logs, `ctest-final.log`, copied `LastTest-final.log` and
`application-{fixed,selected}-final.log`. Earlier successful checks preceded the
final impact validation/context edit and are not counted again. No full matrix,
external solver run, installed-binary audit, distributed execution with the new
weights or experimental mixture/impact validation was performed.
