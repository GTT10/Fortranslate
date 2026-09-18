# Serial EB KH source coupling — 2026-09-18

Implementation `c9cf688`, base `79cb78c`.

`advance_mixture_khrt_shedding` provides a serial population helper using the
existing KH release and four-child construction kernels. It accepts per-parent
gas state/Reynolds number and optional eligibility mask, preserves mixed-liquid
composition, allocates IDs above the shared watermark, and publishes pending
inventories, children, component masses and parent gas-energy exchanges together.
Excluded parents retain their pending state. Capacity and ID overflow reject
without publication.

The allocated EB source invokes this helper after source exchange and before
drift. A parent whose updated trajectory intersects the EB skips KH for that
substep. Children inherit sampled parent pressure for impact dispatch. Their
transverse kinetic-energy contribution is deposited through the parent's gas
stencil. The whole interval retains gas/particle/film/ID rollback. The fixed
population interface cannot enable KH children without an ID argument. The
serial application continues to reject KH child mode; reactive and application
connections are subsequent work.

## Completed checks

GNU Fortran 13.3.0, Open MPI 4.1.6, Debug/Ninja; reused fixed tests-enabled and
selected H2O2 tests-disabled configurations. This executed test is serial.

```sh
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-debug --target test_spray_eb_impact_coupling --parallel 2
env -u LD_LIBRARY_PATH -u LIBRARY_PATH ctest --test-dir build/spray-amr-reactive-debug -R '^spray_eb_impact_coupling$' --output-on-failure
env -u LD_LIBRARY_PATH -u LIBRARY_PATH cmake --build build/spray-amr-reactive-selected-debug --target pelef_selected_spray_tree_runtime --parallel 2
python3 -m py_compile tools/check_project_contract.py
python3 tools/check_project_contract.py --source .
git diff --check
```

CTest **1/1 passed**, test 0.17 s, total 0.24 s. Existing drag, evaporation,
heating, deposition, splash and film checks remain active. Added cases use a
1000 m/s mixed parcel away from the wall over 1e-8 s, with drag, heating and
evaporation disabled. Both impact-enabled and specular paths retain the parent
and create IDs 101–104 above retired watermark 100, with parent ID 1 and motion
during the same source interval. A capacity-one failure preserves gas,
temperature, components, population, ledger and watermark.

A 200 m/s crossing case produces four splash children and one film, with no KH
pending inventory or extra KH IDs. Gas/particle/film totals minus wall exchange
close within relative 3e-12, scale floor 1e-30. Fixed and selected library builds,
provenance, Python syntax, contract and whitespace checks passed. No unexpected
failures or skips. Selected runtime, full matrix, external solver parity,
installed audit and physical breakup validation were not run.

Artifacts: `build/spray-eb-khrt-source-checks/{provenance.log,build.log,ctest.log,build-selected.log}`.
