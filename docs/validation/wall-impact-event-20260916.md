# Transactional wall-impact event — 2026-09-16

Source: parent `d64e725` plus this change. GNU Fortran 13.3 Debug,
`build/port-spray-debug`, MPI disabled.

## Scope and explicit closure

The existing classification/children/film constructors are connected with parent
transitions in one standalone transaction. All six regimes return updated parent,
children, film patch and signed wall impulse/energy. Total mass stays in liquid
states. Rebound and thermal-remnant velocities follow the pinned formulas;
retained parent placement uses a one-diameter outward offset from the supplied
impact location. Remaining flight time and CFD boundary/restart integration are
still caller work.

The reference wet classifier sets splashed mass but leaves deposited mass zero.
The new event closes that incomplete branch by depositing 1-ms, without film
entrainment. This deliberate closure is recorded in the reference manifest and
is not exact wet-event reference parity or physical qualification.

## Focused evidence

Three entries passed, zero final failures/skips, in 0.03 s: wall event, state and
correlations. The event test executes all six regimes, verifies liquid mass and
wall momentum/energy closure, retained-parent outward velocity, restitution,
thermal parent/child split, dry/wet child/deposit split and material identity
through the constructors. Failure during deposition after child generation
preserves every caller output and returns zero wall budget.

The first event fixture used normal incidence, Kv=150 and film ratio .03 for
wet splash. The reference normal child-speed expression gives -0.0069 m/s,
so the constructor correctly rejected it. That condition now has an explicit
rejection/rollback assertion; the successful wet event uses Kv=300 where child
velocity is outward. No bounds or reference constants were loosened. The focused
three-entry run passed after the assertion update.

Build, mechanism provenance, project contract and whitespace checks passed.
No Python source changed; no broad configuration matrix was run.
Artifacts: `build/port-spray-debug/wall-event-build.log` and
`build/port-spray-debug/Testing/Temporary/LastTest.log`. Source identities remain
in `references/spray_wall_impact.json`.

```sh
cmake --build build/port-spray-debug --target test_spray_core --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_wall_(event|state|impact)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 tools/check_project_contract.py --source .
git diff --check
```
