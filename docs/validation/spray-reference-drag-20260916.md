# PelePhysics drag option — 2026-09-16

## Source and implementation

Tested parent: `79796b3378efdc5c8edc95fbbd787b6158308cab` plus this commit's
implementation/tests. Source SHA-256 values and build/test logs are stored in
`build/reference-drag-checks/`.

`references/spray_drag.json` pins PelePhysics
`6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`, `Source/Spray/Drag.H`, with its
SHA-256 and license. The new `drag_model='pelephysics'` routes the piecewise
Cd correlation and ETAB multiplier into existing conservative coupled updates.
Schiller-Naumann remains the default. The model selection is validated and
included in restart identity; existing default contexts remain unchanged.

Cd*Re is evaluated directly, avoiding division by a tiny Reynolds number.
The reference zero-relative-velocity branch returns zero rate. Negative
ETAB drag multipliers reject the frozen update without changing the parcel or
feedback. KH-RT and breakup-disabled parcels do not apply ETAB deformation.

## Focused verification

GNU Fortran 13.3.0, existing Debug spray build, MPI disabled:

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d test_spray_core --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(drag|application_cli)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Two CTest entries passed in 26.59 s, including 23 CLI tests; zero failed or
skipped tests. Build, provenance, Python syntax and project contract passed.
The compiler emitted real-comparison and existing bounds-analysis warnings.

Kernel checks cover Re=0, 0.5, 100 and 2000, the deformation multiplier and
negative-multiplier rollback alongside the existing Schiller-Naumann and
conservative exchange checks. The new CLI ETAB case checks system conservation,
byte-identical gas/parcel CSVs and identical final history after split restart,
changed-model restart rejection and invalid selector rejection.

This ports the force correlation. Frozen exponential integration, gas-film
transport and inventory-based source bookkeeping remain PeleF conventions.
No external solver trajectory parity, full-suite test or physical spray
qualification was performed. Mixture droplets, child dispersion and particle
MPI/AMR/EB remain open.
