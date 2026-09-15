# Stationary film heating and evaporation — 2026-09-16

Source: parent `b85f843` plus this change; GNU Fortran 13.3 Debug,
`build/port-spray-debug`, fixed H2/O2 mechanism, MPI disabled.
Reference `WallFilm.H` identity is pinned in `references/spray_wall_impact.json`.

## Implemented behavior

Film gas heat transfer, wall conduction and single-component Chilton-Colburn
mass transfer are connected to the coupled step and CLI. Bounded explicit
substeps use the existing transport/thermodynamic data. Diameter remains fixed;
height follows mass and current density. Film below the cutoff completely
transfers to gas and is removed. Wall heat enters the surface ledger separately;
gas energy is the discrete integrated film-energy loss plus supplied wall heat.
This adapts reference correlations to PeleF conservation conventions rather than
claiming full reference trajectory/source parity. A single input liquid thermal
conductivity is shared across materials. Lateral flow and moving-wall film
physics remain unimplemented.

## Focused evidence

Three CTest entries passed in 23.30 s, zero failures/skips: new film evolution,
wall coupling and CLI. CLI contains 21 tests (22.60 s CTest wall time). The core
case checks thermal rates, explicit temperature update, mass/footprint invariants,
combined gas/film/wall energy, and invalid-conductivity rollback. The CLI verifies
film mass reduction, gas vapor, exact gas/parcel/film CSV and final history on
restart, changed-conductivity rejection, and film disappearance followed by exact
restart with an empty film array. The disappearance fixture deliberately raises
the numerical height cutoff; it is an edge-case test, not a physical evaporation
time prediction. An earlier three-entry run passed before that edge case was added.

The example completed 10 steps to 1e-6 s. Maximum recorded film mass was
9.999772207233305e-11 kg; final film mass 9.998874012886235e-11 kg and temperature
330.10234042131094 K. Maximum relative system mass/energy errors were
4.0741130598780306e-16 and 8.390522635206275e-16.
Build, provenance, Python syntax, project contract and whitespace passed.
No failed checks or broad matrix; no physical film qualification.

Artifacts: `build/wall-film-evaporation-example/` includes outputs, checkpoint,
run log and `result.json`; `build/port-spray-debug/film-evolution-build.log` and
`build/port-spray-debug/Testing/Temporary/LastTest.log`. CLI fixtures are temporary.

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d test_spray_core test_spray_integration --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(film_evolution|application_cli|integration_wall_coupling)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
mkdir -p build/wall-film-evaporation-example
cd build/wall-film-evaporation-example
../port-spray-debug/pelef_spray_3d ../../cases/spray_3d/wall_film_evaporation.nml
```

Use a fresh output directory for the example.
