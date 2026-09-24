# PeleMP vapor-state cap integration — 2026-09-25

Implementation `19e3a093b4f02a27e0a49db4c6cee821e950e10f`, base `b1b05bc`.
Correction: the EB input mentioned below was the single-level driver. The
separate fixed/selected EB tree frontend still lacked the new field at this
revision; [the follow-up](pelemp-eb-tree-input-20260925.md) adds it and records
changed-rank execution. The historical build successes below did not prove
that new input was accepted by the EB tree frontend.
Frozen source identities and the distinction between Watson boiling latent heat
and ordinary EOS vapor-state latent heat are in
[the source audit](spray-property-scope-20260925.md).

## Change

The opt-in `pelemp-cap` liquid mode stores a separate normal boiling temperature.
Watson supplies the pressure-dependent boiling cap. Clausius saturation uses
gas enthalpy at the droplet temperature minus stored liquid energy at the
capped surface temperature; Antoine also uses that surface temperature.
Constant-cp liquid energy reproduces the pinned PeleMP latent formula;
polynomial cp retains PeleF's existing integrated energy convention.

Single-liquid and mixed-liquid film builders supply the required gas enthalpy.
The single-liquid accepted-state saturation check forwards its existing EOS
enthalpy. Mixture mechanics' temporary effective liquids clear the component
boiling field along with the existing boiling model. Stored liquid energy,
transfer conservation and default models are preserved.

Regular scalar, per-port and mixture inputs expose the new field and bind it
into restart context. MPI AMR and single-level EB mixture namelists also forward and serialize
it. Regular default context is unchanged; expanded AMR/EB namelist serialization
can invalidate older build contexts. No cross-build compatibility claim is made.

## Completed checks

GNU Fortran 13.3.0, Debug/Ninja, Python 3.12.3; existing fixed tests-enabled and
selected H2/O2 tests-disabled configurations reused unchanged. Fixed core,
regular spray, MPI AMR spray and MPI EB spray built. All three corresponding
selected applications built; MPI ABI audit passed. Compiler warnings remain.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_core pelef_spray_3d pelef_mpi_amr_spray_3d pelef_mpi_eb_spray_3d --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_spray_3d_selected pelef_mpi_amr_spray_3d_selected pelef_mpi_reactive_eb_tree_3d_selected_spray --parallel 2
env -u DISPLAY -u LD_LIBRARY_PATH -u LIBRARY_PATH HWLOC_COMPONENTS=-x86 OMPI_MCA_hwloc_base_binding_policy=none ctest --test-dir build/spray-amr-reactive-debug -R '^spray_(pelemp_cap|surface_cap|mixture_rates|evaporation)$' --output-on-failure --output-log build/pelemp-cap-checks/ctest-final.log
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/spray/test_spray_application.py --exe build/spray-amr-reactive-debug/pelef_spray_3d SprayApplication.test_pelemp_cap_restart SprayApplication.test_boiling_surface_cap_restart
env -u LD_LIBRARY_PATH -u LIBRARY_PATH python3 tests/spray/test_spray_application.py --exe build/spray-amr-reactive-selected-debug/pelef_spray_3d_selected SprayApplication.test_pelemp_cap_restart
python3 -m py_compile tests/spray/test_spray_application.py tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Final CTest **4/4 passed**, total **0.06 s**. The new case checks capped Clausius
latent heat, Raoult partitioning, failure without enthalpy, atomic mixture output,
capped Antoine fit bounds and invalid normal boiling temperature. Existing
surface-cap, evaporation and mixture-rate checks passed unchanged.

Fixed application checks **2/2 passed in 2.271 s**; selected application check
**1/1 passed in 1.728 s**. The new method runs scalar, per-port and mixed-liquid
variants in each executable. Each compares uninterrupted and step-three restart
gas/parcel CSV bytes exactly, rejects a changed normal boiling point and an
invalid supercritical boiling point, and checks existing mass/energy/momentum
bounds (relative 2e-12 for mass/energy, absolute 1e-18 for momentum). Fixtures
disable gas chemistry and run 2-by-2-by-2 cells to 1 microsecond. They exercise
synthetic liquid inputs, not calibrated physical properties.

Initial checks failed because the new unit fixture reused a component ID and
the per-port CLI fixture omitted its count and then geometry. Correcting the
fixtures preserved the existing rejection rules and tolerances. The core test
target was rebuilt after the unit correction. An edit helper also exited on a
Python syntax error before writing; the subsequent patch succeeded. Final
provenance, syntax, project-contract and whitespace checks passed, with no skips.

Artifacts: `build/pelemp-cap-checks/`, including fixed/selected build logs,
`build-test-final.log`, `ctest-final.log`, `LastTest-final.log`, and
`application-{fixed,selected}-final.log`. Initial failure logs are retained there.
CLI fixture directories are temporary and removed by the existing test harness.

MPI AMR/EB were compiled but not executed with the new option. No full matrix,
installed-binary audit, external solver run, convergence, experimental boiling
or full spray-model parity is established by these checks.
