# Per-port liquid property models — 2026-09-15

Source: parent `ecce99a` plus this change. GNU Fortran 13.3 Debug,
`build/port-spray-debug`, fixed H2/O2 mechanism, MPI disabled.

## Scope and boundary

Connect per-port polynomial density, bounded Antoine saturation pressure and
Clausius near-boiling vapor-state caps to the existing material pool. Each port
has independent model selectors, coefficients, temperature ranges and critical
temperature. The established kernels and validity bounds are unchanged.
Property settings alone enable the pool. They are validated before injection
and encoded in restart identity; inactive slots reject nondefault values.
Port models default independently of global models to avoid accidental reuse
of another liquid's fit. Mixture droplets and per-material breakup controls
remain outside this increment.

## Focused verification

Four CTest entries passed, zero failures/skips, in 14.06 s: density, Antoine,
surface cap and application CLI. The CLI contains 17 tests and passed in
14.03 s (CTest wall time). The new fixture combines 290 K water using polynomial
density and Antoine pressure with 400 K water using constant density and a
Clausius cap. It checks final mass/density/diameter consistency, conservation,
exact gas/parcel CSV restart and final history equality. Density, Antoine and
critical-temperature changes reject restart. Out-of-range Antoine injection,
an incompatible cap and unused port coefficients also reject.

Mechanism provenance, Python syntax, project contract and whitespace checks
passed. A first build rejected a new line over 132 columns; wrapping it resolved
the compile error. No broad matrix or physical qualification was run.

Artifacts: `build/port-spray-debug/port-fits-build.log` and
`build/port-spray-debug/Testing/Temporary/LastTest.log`. CLI fixture directories
are temporary. Source mechanism SHA-256 remains
`0efc6c52862741a29e0c29b65d979c7d8cb409db5282bca83b9c5437b3d8c8d4`.

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(application_cli|density|antoine|surface_cap)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```
