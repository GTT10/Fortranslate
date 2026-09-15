# Per-port size distributions — 2026-09-16

## Implementation and provenance

Tested parent: `8c5368e2ee6ca17d79a2016b6014925ad2615b7a` plus this commit's
implementation/tests. Exact source hashes and logs are in
`build/port-size-checks/`.

The pinned PelePhysics SprayJet owns its own distribution instance.
`references/spray_size_distribution.json` now records that source alongside
the distribution laws. Fortran `spray_size_options` records let the multi-port
injector select a complete size configuration per port. CLI overrides cover
law, mass/number weighting, standard deviation, shape and diameter bounds;
blank/-1 values inherit global fields independently.

The application resolves and probes all ports before running, passes their
settings into initial and shared/independent scheduled injections, and adds
explicit overrides to restart context. Other liquid properties/material IDs
are independent of this choice. Default configuration retains prior contexts.

## Focused verification

GNU Fortran 13.3.0, existing Debug spray build, MPI disabled:

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d test_spray_core --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(ports|application_cli)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Two CTest entries passed in 34.50 s, including 27 CLI tests, with no failed
or skipped tests. An initial compile failed because the optional options-array
declaration was inserted in the wrong procedure; its placement was corrected
before the successful build. Compiler warnings include real comparisons and
existing bounds-analysis diagnostics. Provenance, Python syntax and project
contract checks passed.

The new case combines a number-weighted normal port with a mass-weighted
Weibull port, each with different width/shape/bounds. Both shared and independent
schedules are checked against Python NormalDist and inverse-Weibull diameter
oracles, their respective equal-count/equal-mass invariants and system
conservation. Split restart gives byte-identical gas/parcel CSVs and identical
final history. Changed width rejects restart; invalid and unused-port settings
are rejected. Existing global distribution and default multi-port tests pass.

## Remaining boundary

This connects independent diameter statistics, not mixture droplets or
physical breakup-child formation/dispersion. Deterministic quadrature and
finite mass normalization remain distinct from reference random sampling.
No full-suite, external solver parity or physical validation was run.
