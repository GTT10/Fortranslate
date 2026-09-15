# Multiple pure-liquid populations — 2026-09-15

Source: parent `f70ebb0` plus this change; GNU Fortran 13.3 Debug,
`build/port-spray-debug`, fixed H2/O2 gas mechanism, MPI disabled.

## Scope

Each port can select its pure-liquid vapor species, density, heat capacity,
reference temperature/pressure and latent heat. Parcel material IDs route
injection mass, liquid energy and evaporated species through both coupling
half steps, shared/independent schedules, output and opt-in checkpoint records.
The application checks restart IDs against the configured material pool before
computing inventories. Legacy single-liquid formats remain unchanged.

Application overrides currently require constant density, Clausius pressure
and rejection of boiling. Breakup viscosity/surface tension remain shared.
These are separate pure-liquid populations, not composition-varying droplets.
H2O2 property values in the fixture are illustrative, with chemistry disabled.

## Focused evidence

Seven CTest entries passed, zero failures/skips, in 17.74 s:
port states, ports, coupled integration, rollback, checkpoint, materials and
application CLI. The CLI entry contains 16 passing tests (12.94 s).
The material integration case checks both evaporated species, their separate
inventories, total mass/momentum/energy, invalid material rollback, legacy schema
rejection and material roundtrip combined with breakup/split/RT/force fields.
The CLI case checks different densities, both schedules, exact gas and parcel
CSV restart and final history row, changed-property rejection, invalid vapor,
invalid density and unused material input rejection.

The supplied case completed 20 steps to 2e-6 s with 40 parcels of materials 1/2.
Injected mass: 3.5999999999999984e-10 kg. Maximum relative mass error:
9.92720400706839e-16; energy error: 2.0581939770169857e-15.
This is conservation/integration evidence, not physical spray qualification.

Mechanism source SHA-256 is
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.
Provenance, Python syntax, project contract and whitespace checks passed.
An initial build failed on a new line exceeding 132 columns; it was wrapped
and the build passed. An example launch from the wrong directory failed before
execution; the corrected command below completed.

Artifacts: `build/port-materials/` contains outputs, checkpoint, `run.log`,
`failed-launch.log` and `result.json`; build output is
`build/port-spray-debug/materials-build.log`. CTest details are in
`build/port-spray-debug/Testing/Temporary/LastTest.log`; CLI cases use temporary
directories.

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d test_spray_integration test_spray_core --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(integration_(materials|coupled|rollback|checkpoint)|ports|port_states|application_cli)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
mkdir -p build/port-materials
cd build/port-materials
../port-spray-debug/pelef_spray_3d ../../cases/spray_3d/port_materials.nml
```

Use a fresh output directory for the example.
