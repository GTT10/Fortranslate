# Number-weighted Weibull injection — 2026-09-16

## Source and change

Tested parent: `3f5b711c8cf90adde98dd62be084ca059fd49f1c` plus this commit's
implementation/tests. Tested file SHA-256 values and build/test logs are in
`build/number-cdf-checks/`.

PelePhysics `6e8026fc6ed121c07e32c9d3fe1d1c919bafc64b`,
`Source/Spray/Distribution/Distributions.cpp`, samples the Weibull inverse CDF
as a droplet-number distribution. `references/spray_size_distribution.json`
records the SHA-256 and license. PeleF previously exposed only a mass-CDF
interpretation with equal parcel mass.

`diameter_weighting='number'` now assigns equal represented droplet counts
within each batch, with unequal parcel masses. The common multiplicity is
normalized to the prescribed dose. The existing default `mass` path retains
its calculation. Optional API arguments propagate through initial, shared-port
and independently scheduled injection. The selector is validated and enters
restart identity. Default restart contexts remain unchanged.

## Focused verification

GNU Fortran 13.3.0, existing Debug spray build, MPI disabled:

```sh
cmake --build build/port-spray-debug --target pelef_spray_3d test_spray_core --parallel 8
ctest --test-dir build/port-spray-debug -R '^spray_(size_distribution|application_cli)$' --output-on-failure
python3 tools/check_mechanism_provenance.py --bundle mechanisms/h2o2_full.json
python3 -m py_compile tests/spray/test_spray_application.py
python3 tools/check_project_contract.py --source .
git diff --check
```

Two CTest entries passed in 28.16 s, including 24 CLI tests; no failures or
skips. Build, provenance, syntax and contract checks passed. Build warnings
include real comparisons and existing bounds-analysis diagnostics.

The kernel checks inverse number-CDF quantiles, equal multiplicities, unequal
single-drop masses and total injected mass. Existing mass-CDF checks also pass.
New CLI cases cover single-port, shared two-port and independently scheduled
two-port batches; each checks conservation, equal batch multiplicities,
byte-identical gas/parcel CSVs and identical final history after split restart.
Changed-weighting restart and invalid-selector cases are rejected.

## Scope

Finite deterministic midpoint quadrature, explicit diameter truncation and
exact-dose normalization differ from upstream random sampling. The Weibull
scale is not generally the arithmetic mean. This completes the selectable
number-basis injection path; it does not implement normal/lognormal size
laws, physical breakup-child dispersion, mixture droplets or particle MPI/AMR/EB.
No full-suite, external trajectory parity or physical spray qualification was run.
